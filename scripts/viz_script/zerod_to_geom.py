import sys
import os
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.solver import Solver0D
from src.misc import m2d, d2m
from src.centerlines import Centerlines

import numpy as np
import shutil

try:
    import sv_rom_extract_results as extract_results
except ImportError:
    print('please run using simvascular --python -- this_script.py')
    exit(1)
    


def add_avg_max_mPAP(centerlines: Centerlines, solver: Solver0D):
    #! SPlit this up into 2 fgunctions more similar to results displayed by sv 0D rom
    
    array_num = centerlines.centerlines.GetPointData().GetNumberOfArrays()
    array_names = [centerlines.centerlines.GetPointData().GetArrayName(i) for i in range(array_num)]
    mpap = {}
    for name in array_names:
        if name.startswith('pressure'):
            mpap[float(name.split('_')[-1])] = centerlines.get_pointdata(name)
    
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([d2m(mpap[tidx]) for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'mPAP')
    
    max_pap = solver.inflow.t[np.where(solver.inflow.Q == solver.inflow.max_inflow)][0] + solver.inflow.tc * (solver.simulation_params['number_of_cardiac_cycles'] - 1)

    centerlines.add_pointdata(d2m(mpap[max_pap]), 'maxPAP_' + str(max_pap))
    
    return centerlines

def only_summary(centerlines: Centerlines):
    ''' Deletes non-summary values  '''
    #! MAKE it so it deletes all other values?
    array_num = centerlines.centerlines.GetPointData().GetNumberOfArrays()
    array_names = [centerlines.centerlines.GetPointData().GetArrayName(i) for i in range(array_num)]
    mpap = {}
    for name in array_names:
        if name.startswith('pressure_') or name.startswith('flow_'):
            centerlines.centerlines.GetPointData().RemoveArray(name)
    return centerlines
            

# TODO: WRITE THESE TO MATCH THE FORMULAS            


            

def main(args):
    
    solver_dir = os.path.dirname(args.solver_file)
    solver_filename = os.path.basename(args.solver_file)
    
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
    
    
    params.output_file_name = os.path.splitext(solver_filename)[0] + '_centerline_results'
    
    # if centerlines not in the same dir as the solver
    if os.path.dirname(args.centerlines_file) != solver_dir:
        shutil.copy(args.centerlines_file, solver_dir )
    
    # process
    post = extract_results.Post(params, None)
    post.process()
    # ! can we skip all the extra processing and just get the array.
    
    centerlines = Centerlines()
    centerlines_file = os.path.join(solver_dir, params.output_file_name + '.vtp')
    centerlines.load_centerlines(centerlines_file)
    centerlines = add_avg_max_mPAP(centerlines, solver)
    if args.summary:
        centerlines = only_summary(centerlines)
        
    centerlines.write_centerlines(centerlines_file)
    
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', dest = 'solver_file', help = 'solver file path')
    parser.add_argument('-c', dest = 'centerlines_file', help = 'centerlines_file')
    parser.add_argument('-s', dest = 'summary', action = 'store_true', default = False, help = 'save summary only')
    args = parser.parse_args()
    
    main(args)