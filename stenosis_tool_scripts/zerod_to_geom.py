
from copyreg import pickle
from functions.utils import *
from functions.project_org import ProjectPath
from functions.solver_io import Solver0D
from functions.centerlines import Centerlines
import os
import numpy as np
import shutil

try:
    import sv_rom_extract_results as extract_results
except ImportError:
    print('please run using simvascular --python -- this_script.py')
    exit(1)
    


def add_avg_max_mPAP(centerlines_file, max_pap):
    centerlines = Centerlines()
    centerlines.load_centerlines(centerlines_file)
    array_num = centerlines.centerlines.GetPointData().GetNumberOfArrays()
    array_names = [centerlines.centerlines.GetPointData().GetArrayName(i) for i in range(array_num)]
    mpap = {}
    for name in array_names:
        if name.startswith('pressure'):
            mpap[float(name.split('_')[-1])] = centerlines.get_pointdata(name)
    
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([mpap[tidx]/1333.22 for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'mPAP')
    
    centerlines.add_pointdata(mpap[max_pap]/1333.22, 'maxPAP_' + str(max_pap))
    
    centerlines.write_centerlines(centerlines_file)


def tool_main(args):
    
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
    for bc in solver.bc:
        if bc['bc_type'] == 'FLOW':
            tc = bc['bc_values']['t']
            
    params.time_range = (tc[-1] * (solver.simulation_params['number_of_cardiac_cycles'] - 1), tc[-1] * solver.simulation_params['number_of_cardiac_cycles'])
    print(params.time_range)
    
    # model input
    params.oned_model = None
    params.centerlines_file = os.path.basename(args.centerlines_file)
    
    
    params.output_file_name = os.path.splitext(solver_filename)[0] + '_centerline_results'
    shutil.copy(args.centerlines_file, solver_dir )
    post = extract_results.Post(params, None)
    print(post.geos['cent'])
    post.process()
    
    results = np.load(os.path.join(solver_dir, os.path.splitext(solver_filename)[0] + '_all_results.npy'), allow_pickle=True).item()
    pressures = results['pressure']['P_BC0_inlet_V0'][-1 * solver.simulation_params['number_of_time_pts_per_cardiac_cycle']:]
    x_max = results['time'][-1 * solver.simulation_params['number_of_time_pts_per_cardiac_cycle']:][np.where(pressures == pressures.max())][0]
    
    
    add_avg_max_mPAP(os.path.join(solver_dir, params.output_file_name + '.vtp'), x_max)
    
    

if __name__ == '__main__':
    
    parser, _, dev, tool = create_parser(description='Converts 0D results to centerlines')
     
   
    tool.add_argument('-solver_file', help = 'solver file path')
    tool.add_argument('-centerlines_file')
     
    args = parser.parse_args()
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        raise NotImplementedError