
from src.data_org import StenosisToolResults, DataPath
from src.misc import create_parser, get_basename
from src.file_io import write_json, check_exists
from src.solver import Solver0D
import numpy as np
import os



#############
# functions #
#############
    
def compute_stenosis_vessels(solver: Solver0D, plausible_gens = {0,1,2,3}, occlusion_range = (.5, .75)):
    ''' Uses a poisson distribution (lambda = len(plausible_gens)) to determine number of stenosis point to actually use'''
    
    # get total number of possible stenosis points for generations 0 - 3
    all_stenosable_vessels = []
    branch_tree = solver.get_branch_tree()
    for node in solver.tree_bfs_iterator(branch_tree):
        if node.generation in plausible_gens:
            all_stenosable_vessels.append(node)

    print('Total number of stenosable branch:', len(all_stenosable_vessels))
    n = min(np.random.poisson(lam = 2 * len(plausible_gens)), len(all_stenosable_vessels) )
    n_vessels = np.random.choice(all_stenosable_vessels, n, replace = False)
    occlusions = np.random.uniform(occlusion_range[0], occlusion_range[1], size = n)
    
    return n_vessels, occlusions


def create_new_solver(old_solver: Solver0D, out_dir, n_vess, occlusions):
    
    vessels = []
    old_sc = []
    old_rp = []
    old_l = []
    old_c = []
    exp_occlusions = []
    
    for idx in range(len(n_vess)):
        for vidx in range(len(n_vess[idx].vess_id)):
            v = n_vess[idx].vessel_info[vidx]
            vessels.append(n_vess[idx].vess_id[vidx])
            old_sc.append(v['zero_d_element_values']['stenosis_coefficient'])
            old_rp.append(v['zero_d_element_values']['R_poiseuille'])
            old_l.append(v['zero_d_element_values']['L'])
            old_c.append(v['zero_d_element_values']['C'])
            exp_occlusions.append(occlusions[idx])
    
    exp_occlusions = np.array(exp_occlusions)
    new_sc = new_sten_coeff(np.array(old_sc), occlusion= exp_occlusions)
    new_rp = new_r_poiseuille(np.array(old_rp), occlusion= exp_occlusions)
    new_l = new_inductance(np.array(old_l), occlusion= exp_occlusions)
    new_c = new_capacitance(np.array(old_c), occlusion= exp_occlusions)
    
    for idx in range(len(n_vess)):
        for v in n_vess[idx].vessel_info:
            v['zero_d_element_values']['stenosis_coefficient'] = new_sc[idx]
            v['zero_d_element_values']['R_poiseuille'] = new_rp[idx]
            v['zero_d_element_values']['L'] = new_l[idx]
            v['zero_d_element_values']['C'] = new_c[idx]
        
    stenosis_file = os.path.join(out_dir, 'stenosis_vessels.dat')
    print('Changed vessels:', vessels)
    changes = {'all_changed_vessels':vessels,
                'occlusions': exp_occlusions.tolist(),
                               'r_poiseuille_old': old_rp,
                               'r_poiseuille_new': new_rp.tolist(),
                               'sten_coeff_old': old_sc,
                               'sten_coeff_new': new_sc.tolist(),
                               'inductance_old': old_l,
                               'inductance_new': new_l.tolist(),
                               'capacitance_old': old_c,
                               'capacitance_new': new_c.tolist()}
    write_json(stenosis_file, changes)
     
    new_solver_file = os.path.join(out_dir, get_basename(old_solver.solver_file) + '_art_sten.in') 
    
    old_solver.write_solver_file(new_solver_file)


def new_inductance(old_ind, occlusion):
    
    return old_ind / (1 - occlusion)

def new_capacitance(old_c, occlusion):
    return (1 - occlusion) * old_c

def new_sten_coeff(old_sten_coeff, occlusion):
    
    remaining_area_pct = 1 - occlusion
    
    # modifications to a_0 and a_s
    return old_sten_coeff / (remaining_area_pct**2)
    

def new_r_poiseuille(old_r_p, occlusion):

    remaining_area_pct = 1 - occlusion
    radius_change = np.sqrt(remaining_area_pct)
    return old_r_p / (radius_change ** 4)
    
    
###########
# Drivers #
###########

def tool_main(args):
    raise NotImplementedError
    # TODO: CREATE A MODEL/RESULTS CORRESPONDING and run for healthy

def dev_main(args):
    
    org = DataPath(args.root)
    
    plausible_gens = set(args.gens)
    occlusion_range = (args.occlusion[0], args.occlusion[1])
    print('Plausible generations:', plausible_gens)
    print(f'Occlusion Range: {occlusion_range[0]} - {occlusion_range[1]}')
    
    for model_name in args.models:
        model = org.find_model(model_name)
        
        if model.type != 'healthy':
            print(model_name + ' is not a healthy model, therefore no artificial stenosis is needed: skipping.')
            continue
        
        model_results = StenosisToolResults(args.root, model)
        
        # load in base solver file
        main_solver = Solver0D()
        main_solver.read_solver_file(model_results.base_solver)
        
        
       
        for i in range(args.n_versions):
            # compute which vessels
            n_vessels, occlusions = compute_stenosis_vessels(solver = main_solver,
                                    plausible_gens=plausible_gens,
                                    occlusion_range=occlusion_range)
            
            version_dir = check_exists(os.path.join(model_results.artificial_sten_dir, str(i)), mkdir = True)
            
            
            create_new_solver(main_solver,version_dir,  n_vessels, occlusions)
if __name__ == '__main__':
    
    parser, dev, tool = create_parser(desc = 'Generates artificial stenosis files')
    
    
    # dev params
    dev.add_argument('-models', dest = 'models', action = 'append', default = [], help = 'Specific models to run')
    dev.add_argument('-n_versions',type = int, default = 1, help = 'number of stenosis versions')
    dev.add_argument('-occlusion', type = float, default = [.75, .9], nargs = 2, help = 'Occlusion range')
    dev.add_argument('-gens', type = int, default = [0, 1, 2, 3], nargs = '*', help = 'Generations to use')
    args = parser.parse_args()
    
    
        
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)