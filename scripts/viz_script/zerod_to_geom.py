import sys
import os
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.solver import Solver0D
from src.misc import m2d, d2m
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
    for name in array_names:
        if name.startswith('pressure'):
            time = round(float(name.split('_')[-1]), 5)
            new_name = 'pressure_' + str(time)
            centerlines.rename_pointdata(name, new_name)
            mpap[time] = centerlines.get_pointdata(new_name)
            # convert to mmHg
            centerlines.add_pointdata(d2m(mpap[time]), new_name)
    
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([d2m(mpap[tidx]) for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'mPAP')
    
    max_pap = round(solver.inflow.t[np.where(solver.inflow.Q == solver.inflow.max_inflow)][0] + solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1), 5)
    min_pap = round(solver.inflow.t[np.where(solver.inflow.Q == solver.inflow.min_inflow)][0] + solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1), 5)

    centerlines.add_pointdata(d2m(mpap[max_pap]), 'maxPAP_' + str(max_pap))
    centerlines.add_pointdata(d2m(mpap[min_pap]), 'minPAP_' + str(min_pap))
    
    return centerlines

def only_summary(centerlines: Centerlines):
    ''' Deletes non-summary values  '''

    array_names = centerlines.get_pointdata_arraynames()
    for name in array_names:
        if name.startswith('pressure_') or name.startswith('flow_'):
            centerlines.centerlines.GetPointData().RemoveArray(name)
    return centerlines  

def mm_to_cm(centerline_file):
    c = Centerlines()
    c.load_centerlines(centerline_file)
    print(centerline_file)
    path = c.get_pointdata(c.PointDataFields.PATH)
    path*= .1
    c.add_pointdata(path, c.PointDataFields.PATH)
    c.write_centerlines(centerline_file)
    
    
            

def main(args):
    
    solver_dir = os.path.dirname(args.solver_file)
    solver_filename = os.path.basename(args.solver_file)
    # copy over centerliines
    shutil.copy(args.centerlines_file, solver_dir )
    centerline_file_name = os.path.basename(args.centerlines_file)
    centerline_file_path = os.path.join(solver_dir, centerline_file_name)

    
    if not args.cm:
        print('Converting to mm...', end = '\t', flush = True)
        mm_to_cm(centerline_file_path)
        print('Done')
    
    solver = Solver0D()
    solver.read_solver_file(args.solver_file)
    
    params = extract_results.Parameters()
    params.data_names = ['flow','pressure']
    params.output_directory = solver_dir
    params.results_directory = solver_dir
    
    ## Solver parameters.
    params.solver_file_name = solver_filename
    params.model_order = 0
            
    params.time_range = (solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1), solver.inflow.tc * solver.simulation_params['number_of_cardiac_cycles'])
    print(params.time_range)
    
    # model input
    params.oned_model = None
    params.centerlines_file = os.path.basename(args.centerlines_file)
    
    output_file_name = os.path.splitext(solver_filename)[0] + '_centerline_results'
    params.output_file_name = output_file_name
    
    # process
    post = extract_results.Post(params, None)
    print('Writing each timestep...', end = '\t', flush = True)
    post.process()
    print('Done')
    # ! can we skip all the extra processing and just get the array.
    
    print('Adding summary...', end = '\t', flush = True)
    centerlines = Centerlines()
    new_centerlines_file = os.path.join(solver_dir, output_file_name + '.vtp')
    centerlines.load_centerlines(new_centerlines_file)
    centerlines = add_summary_PAP(centerlines, solver)
    centerlines.write_centerlines(new_centerlines_file)
    print('Done')
    
    if args.summary:
        print('Writing Summary...', end = '\t', flush = True)
        centerlines = only_summary(centerlines)
        print('Done')
        
    centerlines.write_centerlines(os.path.join(solver_dir, params.output_file_name + '_summary.vtp'))
    
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', dest = 'solver_file', help = 'solver file path')
    parser.add_argument('-c', dest = 'centerlines_file', help = 'centerlines_file')
    parser.add_argument('-s', dest = 'summary', action = 'store_true', default = False, help = 'save summary only file')
    parser.add_argument('-mm', dest = 'cm', action = 'store_false', default = True, help = 'use flag if model is a mm based model.')
    args = parser.parse_args()
    
    main(args)