
from functions.project_org import ProjectPath
from functions.utils import *
from functions.solver_io import Solver0D
from functions.find_generations import *
import numpy as np
import svzerodsolver as zerod
import matplotlib.pyplot as plt
import json



#############
# functions #
#############


def get_result_file(solver_file):
    return os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_all_results.npy'

def validate_rez(solver_file, waveform_name):
    
    solver = Solver0D()
    solver.read_solver_file(solver_file)
    for bc in solver.bc:
        if bc['bc_name'] == 'INFLOW':
            break
    inflow_tc = bc['bc_values']['t'][-1]
    num_pts = int(solver.simulation_params['number_of_time_pts_per_cardiac_cycle'])
        
    
    sim_results = np.load(get_result_file(solver_file), allow_pickle = True).item()
    
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(1,1 ,figsize=(15, 10))

    start = -3 *  num_pts
    ax.plot(sim_results['time'][start:], sim_results['pressure']['P_BC0_inlet_V0'][start:]/1333.22)
    time = sim_results['time'][-1*num_pts:]
    time = time - time[0] + sim_results['time'][0]
    mpap_sim = np.trapz(sim_results['pressure']['P_BC0_inlet_V0'][-1*num_pts:], time) / inflow_tc
    ax.set_title(f"Inlet Pressure: mPAP = {mpap_sim/1333.22}")

    ax.set_xlabel('time (s)')
    ax.set_ylabel('pressure (mmHg)')
    

    fig.savefig(os.path.join(os.path.dirname(solver_file), waveform_name))

    
def compute_stenosis_vessels(generations, plausible_gens = [0,1,2,3], occlusion_range = (.5, .75)):
    ''' Uses a poisson distribution (lambda = len(plausible_gens)) to determine number of stenosis point to actually use'''
    
    # get total number of possible stenosis points for generations 0 - 3
    total_n = sum([len(generations[x]) for x in plausible_gens])
    all_stenosable_vessels = []
    for gens in plausible_gens:
        branches = generations[gens]
        for br, vessels in branches.items():
            vessid = vessels[index_most_stenosed(vessels)]['vessel_id']
            all_stenosable_vessels.append(vessid)
    print('Total number of stenosable vessels:', total_n)
    n = min(np.random.poisson(lam = 2 * len(plausible_gens)), total_n)
    n_vessels = np.random.choice(all_stenosable_vessels, n, replace = False)
    occlusions = np.random.uniform(occlusion_range[0], occlusion_range[1], size = n)
    
    return n_vessels, occlusions


def create_new_solver(old_solver_file, new_solver_file, n_vess, occlusions):
    
    copy_solver = Solver0D()
    copy_solver.read_solver_file(old_solver_file)
    
    
    old_sc = []
    old_rp = []
    old_l = []
    old_c = []
    for idx in range(len(n_vess)):
        v = copy_solver.get_vessel(n_vess[idx])
        old_sc.append(v['zero_d_element_values']['stenosis_coefficient'])
        old_rp.append(v['zero_d_element_values']['R_poiseuille'])
        old_l.append(v['zero_d_element_values']['L'])
        old_c.append(v['zero_d_element_values']['C'])
    
    new_sc = new_sten_coeff(np.array(old_sc), occlusion= occlusions)
    new_rp = new_r_poiseuille(np.array(old_rp), occlusion= occlusions)
    new_l = new_inductance(np.array(old_l), occlusion= occlusions)
    new_c = new_capacitance(np.array(old_c), occlusion= occlusions)
    
    for idx in range(len(n_vess)):
        v = copy_solver.get_vessel(n_vess[idx])
        v['zero_d_element_values']['stenosis_coefficient'] = new_sc[idx]
        v['zero_d_element_values']['R_poiseuille'] = new_rp[idx]
        v['zero_d_element_values']['L'] = new_l[idx]
        v['zero_d_element_values']['C'] = new_c[idx]
        
        
    copy_solver.description = {'all_changed_vessels': n_vess.tolist(),
                               'r_poiseuille_old': old_rp,
                               'r_poiseuille_new': new_rp.tolist(),
                               'sten_coeff_old': old_sc,
                               'sten_coeff_new': new_sc.tolist(),
                               'inductance_old': old_l,
                               'inductance_new': new_l.tolist(),
                               'capacitance_old': old_c,
                               'capacitance_new': new_c.tolist()}
    copy_solver.write_solver_file(new_solver_file)


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

def dev_main(args):
    
    org = ProjectPath(args.root)
    
    for model_name in args.models:
        model = org.find_model(model_name)
        
        # load in solver file
        main_solver = Solver0D()
        main_solver.read_solver_file(model.model_solver)
        
        generations = find_generations(main_solver)
        plausible_gens = [0, 1, ]
        
        if args.sanity_check:
            # Run a simulation to test
            
            if os.path.exists(os.path.join(os.path.dirname(model.model_solver), f'{model_name}_model_waveforms.png')):
                print('Original Model sanity check already performed.')
            else:
                blockPrint() # silence
                zerod.solver.set_up_and_run_0d_simulation(zero_d_solver_input_file_path=model.model_solver,
                                                last_cycle=False,
                                                save_results_all=True,
                                                save_results_branch=False,
                                                use_steady_soltns_as_ics = True)
                enablePrint()
                validate_rez(solver_file=model.model_solver, waveform_name=f'{model_name}_model_waveforms.png')

        for i in range(args.n_versions):
            # compute which vessels
            n_vessels, occlusions = compute_stenosis_vessels(generations=generations,
                                    plausible_gens=plausible_gens,
                                    occlusion_range=(.75, .9))
            
            
            model_dir_name = model.params['metadata']['name'] + '_' + str(i)
            model_dir = os.path.join(model.artificial_stenosis_dir, model_dir_name)
            if not os.path.exists(model_dir):
                os.mkdir(model_dir)
            sten_model_name = os.path.join(model_dir, model_dir_name +  '.in')
            
            print('model_name:', os.path.basename(sten_model_name))
            print('n_vessels:', n_vessels)
            print()
            create_new_solver(model.model_solver, sten_model_name, n_vessels, occlusions)
            
            sten_txtfile = os.path.join(model_dir, model.STENOSIS_FILE )
            with open(sten_txtfile, 'w') as sten_file:
                json.dump({'vessel_ids': n_vessels.tolist(), 'occlusions': occlusions.tolist()}, sten_file,  indent = 4, sort_keys=True)
            
            
            if args.sanity_check:
                # Run a simulation to test
                blockPrint() # silence
                zerod.solver.set_up_and_run_0d_simulation(zero_d_solver_input_file_path=sten_model_name,
                                                last_cycle=False,
                                                save_results_all=True,
                                                save_results_branch=False,
                                                use_steady_soltns_as_ics = True)
                enablePrint()
                validate_rez(solver_file=sten_model_name, waveform_name=f'{os.path.splitext(os.path.basename(sten_model_name))[0]}_waveforms.png')


if __name__ == '__main__':
    
    parser, _, dev, tool = create_parser(description='Generates artificial stenosis files')
    
    
    # dev params
    dev.add_argument('-root', dest = 'root', type = str, default = '.',  help = 'Root to entire project')
    dev.add_argument('-models', dest = 'models', action = 'append', default = [], help = 'Specific models to run')
    dev.add_argument('-n_versions',type = int, default = 10, help = 'number of stenosis versions')
    dev.add_argument('-no_sanity_check',  dest = 'sanity_check', action='store_false', default = True,  help = 'whether to run a 0D simulation for a sanity check')
    args = parser.parse_args()
    
    
        
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)