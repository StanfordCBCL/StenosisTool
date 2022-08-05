# File: correct_healthy_3D_deprecated.py
# File Created: Sunday, 31st July 2022 2:04:51 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 5th August 2022 1:33:34 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Copy over 3D files to create an appropriate directory to run 3D simulations in. Currently deprecated due to unresolved issues
#!      NEED TO RERUN PREFILE w/o DISTAL PRESSURES SINCE IT IS ADDITIVE.... JUST rewrite prefile and dont copy everything over

from src.bc import BoundaryConditions
from src.file_io import parse_face_names
from src.misc import create_tool_parser, get_solver_path
from src.solver import Solver0D
import argparse
import os
import re
import shutil


def read_pre_file(pre_file):
    
    mapping = {}
    with open(pre_file, 'r') as pfile:
        for line in pfile:
            if line.startswith('set_surface_id_vtp'):
                line = line.split()
                face_id = int(line[2])
                face_name = os.path.splitext(os.path.basename(line[1]))[0]
                mapping[face_name] = face_id
                
    return mapping

def correct_solver_file(solver3d_file, rcr_surfaces, num_ts, ts_size, save_occurence):
    with open(solver3d_file, 'r') as sfile:
        file = ''
        for line in sfile:
            
            if re.search('Number of R.* Surfaces:', line):
                file += 'Number of RCR Surfaces: ' + str(len(rcr_surfaces)) + '\n'
            elif re.search('List of R.* Surfaces', line):
                file += 'List of RCR Surfaces: ' + ' '.join([str(x) for x in rcr_surfaces]) + '\n'
            elif re.search('R.* Values', line):
                file += 'RCR Values From File: True\n'
            elif re.search('Number of Timesteps between Restarts:', line):
                file += f'Number of Timesteps between Restarts: {save_occurence}\n'
            elif re.search('Number of Timesteps:', line):
                file += f'Number of Timesteps: {num_ts}\n'
            elif re.search('Time Step Size:', line):
                file += f'Time Step Size: {ts_size}\n'
            else:
                file += line
        
    with open(solver3d_file,'w') as sfile:
        sfile.write(file)

def main(args):
    # copy over
    for solver_dir in args.solver_dirs:
        print(f'Creating 3D files for {solver_dir}...')
        join = lambda fp: os.path.join(solver_dir, fp)
        three_d_dir = join('three_d_dir')
        join_3d = lambda fp: os.path.join(three_d_dir, fp)
        
        # get solver file and model name
        solver_file = get_solver_path(solver_dir)
        solver = Solver0D()
        solver.read_solver_file(solver_file)
        model_name = solver.simulation_params['model_name']
        
        # get simulation path
        model_dir = os.path.dirname(solver_dir)
        simulation_path = os.path.join(model_dir, 'Simulations', model_name)
        if not os.path.exists(simulation_path):
            print(f'{simulation_path} does not exist... Skipping.')
            continue
        
        if not args.force and os.path.exists(three_d_dir):
            print(f'{three_d_dir} already exists... Skipping.')
            continue
        
        # copy over files
        shutil.copytree(simulation_path, three_d_dir, dirs_exist_ok=True)
        
        # get svpre file for mapping
        for file in os.listdir(three_d_dir):
            if os.path.splitext(file)[-1] == '.svpre':
                svpre_file = join_3d(file)
           
        # write new rcr file 
        rcr_file = join('rcrt.dat')
        bcs = BoundaryConditions()
        bcs.read_rcrt_file(rcr_file)
        outlet_mapping = read_pre_file(svpre_file)
        bcs.bc_list = sorted(bcs.bc_list, key = lambda x: outlet_mapping[x['faceID']])
        bcs.write_rcrt_file(path = three_d_dir, three_d = True)
        
        # get outlets
        outlets = parse_face_names(join('outlet_face_names.dat'))
        rcr_surfaces = sorted([outlet_mapping[x] for x in outlets])
        
        # convert solver file to correct parameters

        solver_3d_file = join_3d('solver.inp')
        
        num_ts = solver.simulation_params["number_of_cardiac_cycles"] * solver.simulation_params["number_of_time_pts_per_cardiac_cycle"]
        ts_size = solver.inflow.tc/solver.simulation_params["number_of_time_pts_per_cardiac_cycle"]
        save_occ = 25
        
        correct_solver_file(solver_3d_file, rcr_surfaces, num_ts, ts_size, save_occ)
    
        
        

if __name__ == '__main__':
    
    tool = create_tool_parser(desc = 'corrects a 3D solver with the correct values')
    
    tool.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force copy over and restart process')
    args = tool.parse_args()
    main(args)
    