# Performs comparison of Control and Stenosis generation resistances to determine if a segment contains severe to-be-fixed stenoses

import numpy as np
import os

from src.solver import Solver0D
from src.file_io import write_json, check_exists
from src.misc import create_tool_parser, get_solver_path

########################
# Assistance Functions #
########################

def diameter_formula(order, age):
    ''' from https://journals.physiology.org/doi/full/10.1152/ajpheart.00123.2020 
    We assume order is roughly 16 - gen'''
    a = 1.203e-4
    b = 0.3927
    c = 0.2381
    d = 0.001 + a * order * np.exp(b * order ) * (age ** c)
    return d

def get_control_resistance(vessel_gen, model_age, branch_len, viscosity = .04):

    # we assume vessel_order is approximately 16 - vessel_gen since gen 0 = mpa = order 16 
    vessel_order = 16 - vessel_gen
    rad = diameter_formula(vessel_order, model_age) /2
    control_resistance = 8 * viscosity * branch_len / (np.pi * (rad ** 4))
    return control_resistance

def get_avg_resistance(node: Solver0D.BranchNode):
    res = 0
    branch_len = node.get_branch_len()
    for vess in node.vessel_info:
        res += vess['zero_d_element_values']['R_poiseuille'] * vess['vessel_length'] / branch_len
    return res
        
    

def find_stenosis_vessels(model_solver: Solver0D, age, r_threshold, plausible_gens = {0,1,2,3}):
    ''' finds stenosis vessels in test_gens
    r_threshold is the percent higher the test resistance is than the control,
    This should only hold true for the first 4 generations = orders 16 - 12'''
    total_len = 0
    sten_len = 0
    stenosis_vessels = []
    control_vessel_radii = []
    vessel_side = []
    branch_tree = model_solver.get_branch_tree()
    viscosity = model_solver.simulation_params['viscosity']
    
    for node in model_solver.tree_bfs_iterator(branch_tree):
        
        branch_len = node.get_branch_len()
        total_len += branch_len
        gen = node.generation
        if gen in plausible_gens:
            
            ## first check the entire branch avg 
            control_res = get_control_resistance(gen, age,branch_len, viscosity)
            avg_res = get_avg_resistance(node)
            
            # if it is a stenosis point, add the entire branch
            if avg_res > control_res * r_threshold:
                stenosis_vessels += node.vess_id
                vessel_side += [node.side for i in range(len(node.vess_id))]
                control_vessel_radii += [diameter_formula(16 - gen, age = age)/2 for i in range(len(node.vess_id))]
                sten_len += branch_len
            # if whole branch averaged is not a stenosis point, check individual
            else: 
                for vidx in range(len(node.vess_id)):
                    vess = node.vessel_info[vidx]
                    control_res = get_control_resistance(gen, age, vess['vessel_length'], viscosity)
                    res = vess['zero_d_element_values']['R_poiseuille']
                    if res > control_res * r_threshold:
                        sten_len += vess['vessel_length']
                        stenosis_vessels.append(node.vess_id[vidx])
                        vessel_side.append(node.side)
                        control_vessel_radii.append(diameter_formula(16 - gen, age = age)/2)
            
    return stenosis_vessels, control_vessel_radii, total_len, sten_len


################
# Fixed Solver #
################

def compute_radii(viscosity, length, rp):
    r = ((8 * viscosity * length) / (rp * np.pi)) ** (1/4)
    return r

def new_inductance(old_ind, rad_rat):
    return old_ind / (rad_rat ** 2)

def new_capacitance(old_c, rad_rat):
    return rad_rat**2 * old_c

def new_sten_coeff(old_sten_coeff, rad_rat):
    # modifications to a_0 and a_s
    return old_sten_coeff / (rad_rat**4)
    
def new_r_poiseuille(old_r_p, rad_rat):
    return old_r_p / (rad_rat ** 4)

def construct_fixed_solver(sten_model: Solver0D, stenosis_dict, out_dir):
    
    viscosity = sten_model.simulation_params['viscosity']
    for idx in range(len(stenosis_dict['stenosis_vessel_ids'])):
        vid = stenosis_dict['stenosis_vessel_ids'][idx]
        vess = sten_model.get_vessel(vid)
        old_r = compute_radii(viscosity, vess['vessel_length'], vess['zero_d_element_values']['R_poiseuille'])
        new_r = stenosis_dict['control_vessel_radii'][idx]
        r_rat = new_r/old_r
        ele = vess['zero_d_element_values']
        
        
        ele['C'] = new_capacitance(ele['C'],r_rat)
        ele["stenosis_coefficient"] = new_sten_coeff(ele["stenosis_coefficient"], r_rat)
        ele['R_poiseuille'] = new_r_poiseuille(ele['R_poiseuille'], r_rat)
        ele['L'] = new_inductance(ele['L'], r_rat)
        
    out_file = os.path.join(out_dir, os.path.splitext(os.path.basename(sten_model.solver_file))[0]) + '_fixed_stenosis.in'
    sten_model.write_solver_file(out_file)
    

########
# Main #
########


def main(args):
    
    for solver_dir in args.solver_dirs:
        base_solver_file = get_solver_path(solver_dir)
        print(f'Identifying stenosis for {base_solver_file}...', end = '\t', flush = True)

        # get test generations
        sten_solver = Solver0D()
        sten_solver.read_solver_file(base_solver_file)
        
        ## retrieve test model
        model_age = sten_solver.simulation_params['age']
    
        r_threshold = args.r_threshold # 4x  higher
        plausible_gens  = set(args.gens) # first 4 generations
        vessel_ids, control_vessel_radii, total_len, sten_len = find_stenosis_vessels( sten_solver, model_age, r_threshold, plausible_gens) 
            
        stenosis_dict = {
                        'r_threshold': r_threshold ,
                        'stenosis_vessel_ids': vessel_ids,
                        'control_vessel_radii': control_vessel_radii,
                        'total_tree_len': total_len,
                        'total_sten_length': sten_len}
        
        
        write_json(os.path.join(solver_dir, 'stenosis.txt'), stenosis_dict)
        
        # construct a fixed version
        fixed_stenosis_dir = os.path.join(solver_dir, 'fixed_stenosis')
        if args.force or not os.path.exists(fixed_stenosis_dir):
            if not os.path.exists(fixed_stenosis_dir):
                os.mkdir(fixed_stenosis_dir)
            construct_fixed_solver(sten_solver, stenosis_dict, fixed_stenosis_dir)
        print('Done')
            
            
            

if __name__ == '__main__':
    
    parser = create_tool_parser(desc='Generates artificial stenosis files')
    
    # dev params
    parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force a new iteration over the old')
    parser.add_argument('-r_threshold', type = int, default = 4,  help = 'how many x control resistance the stenosis resistance must be to be considered a fixable location')
    parser.add_argument('-gens', type = int, default = [0,1,2,3], nargs = '*', help = 'generations to identify fixable stenosis in')
    args = parser.parse_args()
    
    
        
    main(args)