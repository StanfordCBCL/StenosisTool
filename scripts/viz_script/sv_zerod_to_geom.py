# File: sv_zerod_to_geom.py
# File Created: Wednesday, 27th July 2022 11:28:30 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 5th August 2022 1:36:25 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: For a particular directory (and potentially its children directory) if a centerlines exists and a branch_result.npy file exists, then map the 0D results to centerlines.


import sys
import os
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.solver import Solver0D
from src.misc import create_tool_parser, get_basename, get_solver_path, m2d, d2m, get_solver_name
from src.centerlines import Centerlines

import numpy as np
import shutil
np.testing.suppress_warnings(forwarding_rule='always')

try:
    import sv_rom_extract_results as extract_results
except ImportError:
    print('please run using simvascular --python -- this_script.py')
    exit(1)
    
def add_summary_PAP(centerlines: Centerlines, solver: Solver0D):
    ''' adds max, min, and avg'''
    
    array_names = centerlines.get_pointdata_arraynames()
    mpap = {}
    mq = {}
    for name in array_names:
        if name.startswith('pressure'):
            time = round(float(name.split('_')[-1]), 5)
            new_name = 'pressure_' + str(time)
            centerlines.rename_pointdata(name, new_name)
            mpap[time] = centerlines.get_pointdata(new_name)
            # convert to mmHg
            centerlines.add_pointdata(d2m(mpap[time]), new_name)
        if name.startswith('flow'):
            time = round(float(name.split('_')[-1]), 5)
            new_name = 'flow_' + str(time)
            centerlines.rename_pointdata(name, new_name)
            mq[time] = centerlines.get_pointdata(new_name)
    
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([d2m(mpap[tidx]) for tidx in t])
    Q = np.array([mq[tidx] for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    mQ = np.trapz(Q, t, axis =0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'meanPAP')
    centerlines.add_pointdata(mQ, 'meanQ')
    
    max_pap = round(solver.inflow.t[np.where(solver.inflow.Q == solver.inflow.max_inflow)][0] + solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1), 5)
    min_pap = round(solver.inflow.t[np.where(solver.inflow.Q == solver.inflow.min_inflow)][0] + solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1), 5)

    centerlines.add_pointdata(d2m(mpap[max_pap]), 'sysPAP_' + str(max_pap))
    centerlines.add_pointdata(d2m(mpap[min_pap]), 'diaPAP_' + str(min_pap))
    centerlines.add_pointdata(mq[max_pap], 'sysQ_' + str(max_pap))
    centerlines.add_pointdata(mq[min_pap], 'diaQ_' + str(min_pap))
    
    return centerlines

def only_summary(centerlines: Centerlines):
    ''' Deletes non-summary values  '''

    array_names = centerlines.get_pointdata_arraynames()
    for name in array_names:
        if name.startswith('pressure_') or name.startswith('flow_'):
            centerlines.remove_pointdata(name)
    return centerlines  


def mm_to_cm(centerline_file):
    c = Centerlines()
    c.load_centerlines(centerline_file)
    path = c.get_pointdata(c.PointDataFields.PATH)
    path*= .1
    c.add_pointdata(path, c.PointDataFields.PATH)
    c.write_centerlines(centerline_file)
    
def create_resistance_map(centerlines: Centerlines, solver: Solver0D):
    
    branch_ids = centerlines.get_pointdata(centerlines.PointDataFields.BRANCHID)
    resistances = np.zeros(len(branch_ids))
    
    branch_tree = solver.get_branch_tree()
    #! A very rough estimate by adding resistances. Too lazy to compute each vessel along
    for node in solver.tree_bfs_iterator(branch_tree):
        r_p = 0
        for vess in node.vessel_info:
            r_p += vess['zero_d_element_values']['R_poiseuille']
        resistances[np.where(branch_ids == node.branch_id)] = r_p
    
    centerlines.add_pointdata(resistances, 'Resistances')
            
    return centerlines
            
            
    
    
            

def main(args):
    for solver_dir in args.solver_dirs:
        
        data_set = []
        
        for tup in os.walk(solver_dir):
            sdir = tup[0]
            # skip tuning dir
            if os.path.basename(sdir) == 'tuning_dir':
                continue
            
            # get solver file
            solver_file = get_solver_path(sdir)
            if solver_file is None:
                continue
                
            # get results
            result_file = os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_branch_results.npy'
            if not os.path.exists(result_file):
                print('Simulation for ' + result_file +  ' does not exist.')
                continue
            
            # centerlines
            centerline_file = os.path.join(sdir, 'model_centerlines.vtp')
            if not os.path.exists(centerline_file):
                print('Centerline file ' + centerline_file + ' does not exist')
                
            data_set.append((sdir, solver_file, result_file, centerline_file))
            
            if not args.recursive:
                break
                    
        for (base_dir, solver_file, solver_results_file, centerline_file) in data_set:
            
            print('Mapping ' + solver_results_file + ' to centerlines...')
            solver = Solver0D()
            solver.read_solver_file(solver_file)
            
            # convert centerlines paths to cm to adjust for solver
            units = solver.simulation_params['units']
            if units == 'mm':
                print('Converting centerline paths to cm...', end = '\t', flush = True)
                mm_to_cm(centerline_file)
                print('Done')
            
            
            
            # parameters
            params = extract_results.Parameters()
            params.data_names = ['flow','pressure']
            params.output_directory = base_dir
            params.results_directory = base_dir
    
            params.solver_file_name = os.path.basename(solver_file)
            params.model_order = 0
            
            params.time_range = (solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1), solver.inflow.tc * solver.simulation_params['number_of_cardiac_cycles'])

    
            # model input
            params.oned_model = None
            params.centerlines_file = os.path.basename(centerline_file)
    
            output_file_name = get_basename(solver_file) + '_centerline_results'
            params.output_file_name = output_file_name
    
            # process
            post = extract_results.Post(params, None)
            print('Writing each timestep...', end = '\t', flush = True)
            post.process()
            print('Done')
            
    
            print('Adding summary...', end = '\t', flush = True)
            centerlines = Centerlines()
            new_centerlines_file = os.path.join(base_dir, output_file_name + '.vtp')
            centerlines.load_centerlines(new_centerlines_file)
            centerlines = create_resistance_map(centerlines, solver)
            centerlines = add_summary_PAP(centerlines, solver)
            centerlines.write_centerlines(new_centerlines_file)
            print('Done')
        
        
            if args.summary:
                print('Writing Summary...', end = '\t', flush = True)
                centerlines = only_summary(centerlines)
                summary_centerlines = os.path.join(base_dir, params.output_file_name + '_summary.vtp')
                centerlines.write_centerlines(summary_centerlines)
                print('Done')
    
    
    
    

if __name__ == '__main__':
    
    parser = create_tool_parser(desc = 'Converts 0D to 1D centerlines')
    parser.add_argument('-s', dest = 'summary', action = 'store_true', default = False, help = 'save summary only file')
    parser.add_argument('-r', dest = 'recursive', action = 'store_true', default = False, help = 'run recursively for all subdirs')
    args = parser.parse_args()
    
    main(args)