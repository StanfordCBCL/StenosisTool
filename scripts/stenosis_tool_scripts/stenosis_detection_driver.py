# File: stenosis_detection_driver.py
# File Created: Wednesday, 29th June 2022 11:26:57 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 17th October 2022 6:08:29 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Given a directory of a stenosis model, performs comparison of Control and Stenosis generation resistances to determine if a segment contains severe to-be-fixed stenoses

from asyncore import write
import numpy as np
import os
from collections import defaultdict

from src.lpn import Solver0D
from src.file_io import copy_rel_files, write_json, check_exists
from src.misc import create_tool_parser, get_solver_path
from src.stenosis import new_capacitance, new_inductance, new_r_poiseuille, new_sten_coeff

########################
# Non-Config Functions #
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

def get_true_resistance(node: Solver0D.BranchNode):
    res = 0
    for vess in node.vessel_info:
        res += vess['zero_d_element_values']['R_poiseuille'] 
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
            true_res = get_true_resistance(node)
            
            # if it is a stenosis point, add the entire branch
            if true_res > control_res * r_threshold:
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

####################
# Config Functions #
####################

def gen_config(model_solver: Solver0D, age, r_threshold):
    def d():
        return {'occlusions': [], 'stenosis_length': 0, 'total_length':0}
    generations = defaultdict(d)
    config_file = {}
    config_file['generation_info'] = defaultdict()
    config_file['metadata'] = {'r_threshold': r_threshold, 'model_name': model_solver.simulation_params['model_name']}
    

    branch_tree = model_solver.get_branch_tree()
    viscosity = model_solver.simulation_params['viscosity']
    
    for node in model_solver.tree_bfs_iterator(branch_tree):
        gen = node.generation
        branch_len = node.get_branch_len()
        generations[node.generation]['total_length'] += branch_len
        
        ## first check the entire branch avg 
        control_res = get_control_resistance(gen, age,branch_len, viscosity)
        true_res = get_true_resistance(node)
        
        
        # if it is a stenosis point, add the entire branch
        if true_res > control_res * r_threshold:
            generations[gen]['stenosis_length'] += branch_len
            occlusion = 1 - ((control_res/true_res) ** (1/2))
            generations[gen]['occlusions'].append(occlusion)
        # if whole branch averaged is not a stenosis point, check individual
        else: 
            for vidx in range(len(node.vess_id)):
                vess = node.vessel_info[vidx]
                control_res = get_control_resistance(gen, age, vess['vessel_length'], viscosity)
                res = vess['zero_d_element_values']['R_poiseuille']
                if res > control_res * r_threshold:
                    generations[gen]['stenosis_length'] += vess['vessel_length']
                    occlusion = 1 - ((control_res/res) ** (1/2))
                    generations[gen]['occlusions'].append(occlusion)
    
    for gen, vals in generations.items():
        occlusions = np.array(vals['occlusions'])
        if len(occlusions) == 0:
            occlusions = np.zeros(1)
        config_file['generation_info'][gen] = {'stenosis_length_ratio': vals['stenosis_length']/ vals['total_length'], 'average_occlusion': occlusions.mean(), 'std_occlusion': occlusions.std()}
    
    return config_file

################
# Fixed Solver #
################

def compute_radii(viscosity, length, rp):
    r = ((8 * viscosity * length) / (rp * np.pi)) ** (1/4)
    return r


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
        
        if not args.config:
            # NON-config mode
            
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
                copy_rel_files(solver_dir, fixed_stenosis_dir, exclude_solver = True)
                construct_fixed_solver(sten_solver, stenosis_dict, fixed_stenosis_dir)
            else:
                print('Skipping fixed solver...')
        
        else:
            # Config Mode
            
            # retrieve max number of gens
            branch_tree = sten_solver.get_branch_tree()
            max_gen = 0
            for node in sten_solver.tree_bfs_iterator(branch_tree):
                if node.generation > max_gen:
                    max_gen = node.generation
            plausible_gens = set(list(range(max_gen + 1)))
            # retrieve all stenosable vessels
            config_file = gen_config(sten_solver, model_age, r_threshold = r_threshold)
            
            write_json(os.path.join(solver_dir, args.config_file), config_file)
            
        print('Done')

if __name__ == '__main__':
    
    parser = create_tool_parser(desc='Generates artificial stenosis files')
    
    # dev params
    parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force a new iteration over the old')
    parser.add_argument('-r_threshold', type = float, default = 4,  help = 'how many x control resistance the stenosis resistance must be to be considered a fixable location')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-gens', type = int, default = [0,1,2,3], nargs = '*', help = 'generations to identify fixable stenosis in')
    
    group.add_argument('-config', action = 'store_true', default = False, help = 'Generate a config file for artificial stenosis generation')
    
    parser.add_argument('-config_file', default = 'artificial_stenosis.cfg', help = 'File name to save config to if -config is specified')
    
    args = parser.parse_args()
    
    
        
    main(args)