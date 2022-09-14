# File: artificial_stenosis_driver.py
# File Created: Tuesday, 5th July 2022 3:31:59 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 13th September 2022 10:18:40 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Given a particular directory containing a solver file and a config file containing occlusions mapped from another stenosis, construct a healthy model with an appropriate amount of stenosis




from src.misc import create_tool_parser, get_basename, get_solver_path
from src.file_io import copy_rel_files, write_json, check_exists, read_json
from src.solver import Solver0D
import numpy as np
import os
from collections import defaultdict
import shutil



#############
# functions #
#############


# ONLY WORK DOWNWARDS
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
    
def compute_stenosis_vessels(solver: Solver0D, plausible_gens = {0,1,2,3}, occlusion_range = (.5, .75)):
    ''' Uses a poisson distribution (lambda = len(plausible_gens)) to determine number of stenosis point to actually use'''
    
    # get total number of possible stenosis points for generations 0 - 3
    all_stenosable_vessels = []
    branch_tree = solver.get_branch_tree()
    for node in solver.tree_bfs_iterator(branch_tree):
        if node.generation in plausible_gens:
            all_stenosable_vessels.append(node)

    print('Total number of stenosable branch:', len(all_stenosable_vessels))
    
    
    if plausible_gens == {0, 1}:
        n = len(all_stenosable_vessels)
    else:
        n = min(np.random.poisson(lam = 2 * len(plausible_gens)), len(all_stenosable_vessels) )
    n_vessels = np.random.choice(all_stenosable_vessels, n, replace = False)
    
    lengths = {'total_stenosed_len': 0, 'total_stenosable_len': 0, 'total_tree_len':0}
    for node in n_vessels:
        lengths['total_stenosed_len'] += node.get_branch_len()
    for node in all_stenosable_vessels:
        lengths['total_stenosable_len'] += node.get_branch_len()
    for node in solver.tree_bfs_iterator(branch_tree):
        lengths['total_tree_len'] += node.get_branch_len()
    
        
    print('Lengths:',lengths)

    occlusions = np.random.uniform(occlusion_range[0], occlusion_range[1], size = n)
    
    return n_vessels, occlusions, lengths


def create_new_solver(old_solver: Solver0D, new_solver_file, out_dir, n_vess, occlusions, lengths):
    
    
    new_solver = Solver0D()
    new_solver.read_solver_file(old_solver.solver_file )
    vessels = []
    old_sc = []
    old_rp = []
    old_l = []
    old_c = []
    side = []
    exp_occlusions = []
    
    for idx in range(len(n_vess)):
        branch = n_vess[idx].vessel_info
        side.append(n_vess[idx].side)
        vessels.append(n_vess[idx].vess_id)
        old_sc.append([v['zero_d_element_values']['stenosis_coefficient'] for v in branch])
        old_rp.append([v['zero_d_element_values']['R_poiseuille'] for v in branch])
        old_l.append([v['zero_d_element_values']['L'] for v in branch])
        old_c.append([v['zero_d_element_values']['C'] for v in branch])
        exp_occlusions.append(occlusions[idx])
    
    exp_occlusions = np.array(exp_occlusions)
    new_sc = [list(new_sten_coeff(np.array(old_sc[i]), occlusion= exp_occlusions[i])) for i in range(len(old_sc))]
    new_rp = [list(new_r_poiseuille(np.array(old_rp[i]), occlusion= exp_occlusions[i])) for i in range(len(old_rp))]
    new_l = [list(new_inductance(np.array(old_l[i]), occlusion= exp_occlusions[i])) for i in range(len(old_l))]
    new_c = [list(new_capacitance(np.array(old_c[i]), occlusion= exp_occlusions[i])) for i in range(len(old_c))]
    
    for idx in range(len(n_vess)):
        for v in n_vess[idx].vessel_info:
            v['zero_d_element_values']['stenosis_coefficient'] = new_sc[idx]
            v['zero_d_element_values']['R_poiseuille'] = new_rp[idx]
            v['zero_d_element_values']['L'] = new_l[idx]
            v['zero_d_element_values']['C'] = new_c[idx]
        
    stenosis_file = os.path.join(out_dir, 'stenosis_vessels.dat')
    print('Changed vessels:', vessels)
    changes = {'all_changed_vessels':vessels,
                'lpa/rpa': side,
                'occlusions': exp_occlusions.tolist(),
                               'r_poiseuille_old': old_rp,
                               'r_poiseuille_new': new_rp,
                               'sten_coeff_old': old_sc,
                               'sten_coeff_new': new_sc,
                               'inductance_old': old_l,
                               'inductance_new': new_l,
                               'capacitance_old': old_c,
                               'capacitance_new': new_c,
                               'lengths': lengths}
    write_json(stenosis_file, changes)
    
    old_solver.write_solver_file(new_solver_file)


###################
# Config Function #
###################

def compute_from_config(solver: Solver0D, cfg):
    ''' Uses a normal distribution at each generation to determine number of stenosis point to actually use'''
    
    generations = defaultdict(list)
    branch_tree = solver.get_branch_tree()
    for node in solver.tree_bfs_iterator(branch_tree):
        generations[node.generation].append(node)
        
    def d():
        return {'total_length': 0}
    geninfo = defaultdict(d)
        
    for gen, nodelist in generations.items():
        for node in nodelist:
            geninfo[gen]['total_length'] += node.get_branch_len()
    
    cfg = cfg['generation_info']
    
    # if generation # is greater in current model than cfg, just repeat the last one
    if len(generations.keys()) > len(cfg.keys()):
        last_gen = len(cfg.keys()) - 1
        for i in range(last_gen + 1, last_gen + len(generations.keys()) - len(cfg.keys()) + 1):
            cfg[str(i)] = cfg[str(last_gen)]
    
    # get the vessels
    n_vessels = []
    occlusions = []
    lengths = defaultdict(dict)
    for gen, nodelist in generations.items():
        # shuffle nodelist & compute the branches to be stenosed
        target_length = cfg[str(gen)]['stenosis_length_ratio'] * geninfo[gen]['total_length']
        cur_length = 0
        np.random.shuffle(nodelist)
        idx = 0
        while cur_length < target_length:
            node = nodelist[idx]
            cur_length += node.get_branch_len()
            n_vessels.append(nodelist[idx])
            idx+=1
            # compute the occlusion level, Set a max of 95%, since a 100% occlusion is unreasonable and a min of 0
            occlusion = max(min(np.random.normal(loc = cfg[str(gen)]['average_occlusion'], scale = cfg[str(gen)]['std_occlusion']), .95), 0)
            occlusions.append(occlusion)
        

    
        # compute lengths
        lengths[gen]['total_stenosed_length'] = cur_length
        lengths[gen]['total_stenosable_len'] = geninfo[gen]['total_length']
        
            

    return n_vessels, np.array(occlusions), lengths
   
   
    
###########
# Drivers #
###########

def main(args):
    

    plausible_gens = set(args.gens)
    occlusion_range = (args.occlusion[0], args.occlusion[1])
    
    for solver_dir in args.solver_dirs:
        
        base_solver_file = get_solver_path(solver_dir)
        
        # create the proximal test
        artificial_sten_dir = check_exists(os.path.join(solver_dir, 'artificial_stenosis'), mkdir = True)
        join_art = lambda file: os.path.join(artificial_sten_dir, file)
        version_dir = join_art('proximal')
        proximal_solver_file = os.path.join(version_dir, get_basename(base_solver_file) + '_proximal_sten.in') 
        
        if args.force or not os.path.exists(version_dir):
            # load in base solver file
            proximal_solver = Solver0D()
            proximal_solver.read_solver_file(base_solver_file)
            
            n_vessels, occlusions, lengths= compute_stenosis_vessels(solver = proximal_solver,
                                        plausible_gens={0,1},
                                        occlusion_range=(.75, .75))
            if not os.path.exists(version_dir):
                os.mkdir(version_dir)
            # only if one does not already exist
            copy_rel_files(solver_dir, version_dir, exclude_solver = True)
            create_new_solver(proximal_solver, proximal_solver_file,  version_dir,  n_vessels, occlusions, lengths)
        else:
            print('Skipping proximal solver file writing: proximal solver already exists')
        
        if not args.config:
            # Default version...
            for i in range(args.n_versions):
                # reload in base solver file
                tmp_solver = Solver0D()
                tmp_solver.read_solver_file(base_solver_file)
                
                # compute which vessels
                n_vessels, occlusions, lengths= compute_stenosis_vessels(solver = tmp_solver,
                                        plausible_gens=plausible_gens,
                                        occlusion_range=occlusion_range)
                
                version_dir = join_art(str(i))
                if args.force or not os.path.exists(version_dir):
                    if not os.path.exists(version_dir):
                        os.mkdir(version_dir)
                    copy_rel_files(solver_dir, version_dir, exclude_solver = True)
                    new_solver_file = os.path.join(version_dir, get_basename(base_solver_file) + '_art_sten.in') 
                    create_new_solver(tmp_solver, new_solver_file, version_dir,  n_vessels, occlusions, lengths)
                else:
                    print(f'{version_dir} already exists: skipping')
    
        else:
            # Takes in a config
            cfg = read_json(args.config)
            # Default version...
            for i in range(args.n_versions):
                # reload in base solver file
                tmp_solver = Solver0D()
                tmp_solver.read_solver_file(base_solver_file)
                
                # compute which vessels
                n_vessels, occlusions, lengths= compute_from_config(solver = tmp_solver, cfg = cfg)
                
                
                version_dir = join_art(cfg['metadata']['model_name'] + '_' + str(i))
                if args.force or not os.path.exists(version_dir):
                    print(f'Generating version {version_dir}')
                    if not os.path.exists(version_dir):
                        os.mkdir(version_dir)
                    copy_rel_files(solver_dir, version_dir, exclude_solver = True)
                    new_solver_file = os.path.join(version_dir, get_basename(base_solver_file) + '_art_sten.in') 
                    create_new_solver(tmp_solver, new_solver_file, version_dir,  n_vessels, occlusions, lengths)
                    shutil.copy(args.config, version_dir)
                else:
                    print(f'{version_dir} already exists: skipping')
            

if __name__ == '__main__':
    
    parser = create_tool_parser(desc = 'Generates artificial stenosis files')
    
    # dev params
    parser.add_argument('-n_versions',type = int, default = 1, help = 'number of stenosis versions')
    
    mutex_group = parser.add_mutually_exclusive_group()
    group = mutex_group.add_argument_group()
    
    group.add_argument('-occlusion', type = float, default = [.75, .9], nargs = 2, help = 'Occlusion range')
    group.add_argument('-gens', type = int, default = [0, 1, 2, 3], nargs = '*', help = 'Generations to use')
    mutex_group.add_argument('-config', type = str, help = 'Specify a config file to take in.')
    parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force saving a new proximal version')

    args = parser.parse_args()
    
    main(args)