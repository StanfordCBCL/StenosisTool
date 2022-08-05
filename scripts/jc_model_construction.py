# File: jc_model_construction.py
# File Created: Monday, 18th July 2022 3:40:19 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 5th August 2022 1:38:11 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Constructs a new directory from an base one where everything is identical except the solver file is processed to includes junction coefficients


from src.solver import Solver0D
from src.misc import create_tool_parser, get_basename, get_solver_path
import os
import shutil
from src.file_io import copy_rel_files

def convert_to_jc(solver: Solver0D):
    
    for junc in solver.junctions:
        if junc['junction_type'] != 'internal_junction':
            
            inlet = junc['inlet_vessels']
            outlets = junc['outlet_vessels']
            S0 = junc['areas'][0]
            s_outlets = junc['areas'][1:]
            
            assert len(outlets) == len(s_outlets), 'Number of areas does not match number of outlets'
            
            outlet_areas = list(zip(outlets, s_outlets))
            
            density = solver.simulation_params['density']
            for out, S1 in outlet_areas:
                out_vess = solver.get_vessel(out)
                
                # update with stenosis coefficient of junction
                out_vess['junction_coefficient'] = 1.52 * density * ( (S0/S1 - 1) **2) / (2 * S0**2)
                out_vess['zero_d_element_values']['stenosis_coefficient'] += out_vess['junction_coefficient']

def main(args):
    
    for solver_dir in args.solver_dirs:
        print(f'Constructing JC solver for {solver_dir}...',end = '\t', flush = True)
        # construct a copy of the base directory
        model_dir = os.path.dirname(solver_dir)
        new_jc_dir = os.path.join(model_dir, args.outdir)
        if not os.path.exists(new_jc_dir):
            os.mkdir(new_jc_dir)
        elif not args.force:
            print(new_jc_dir + ' already exists. Skipping')
            continue
        
        # copy files
        copy_rel_files(solver_dir, new_jc_dir, exclude_solver=True)
        
        # rename solver file
        solver_file=get_solver_path(solver_dir)
        solver = Solver0D()
        solver.read_solver_file(solver_file = solver_file)
        
        convert_to_jc(solver)

        new_jc_solver = os.path.join(new_jc_dir, get_basename(solver_file) + '_jc.in')
        solver.write_solver_file(new_jc_solver)
        print('Done')
        
        

if __name__ == '__main__':
    
    tool = create_tool_parser(desc = 'Implementing a new method of stenosis junction handling')
    
    tool.add_argument('-o', dest = 'outdir', default = 'jc_solver_dir', help = 'Dirname of out file')
    tool.add_argument('-f', dest = 'force', default = False, action = 'store_true', help = 'Force override files in existing dir')
    args = tool.parse_args()

    main(args)
    
    