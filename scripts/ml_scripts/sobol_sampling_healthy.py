# File: sobol_sampling_healthy.py
# File Created: Friday, 19th August 2022 4:22:32 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 15th September 2022 8:19:52 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Use Sobol sampling to retrieve a distribution of potential diameter changes for a particular healthy model. Takes in an artificial stenosis directory.



from math import ceil
from unittest import result
import tqdm
import ray
import os
import pandas as pd
import numpy as np
from scipy.stats.qmc import Sobol, scale
import copy
from svzerodsolver.runnercpp import run_from_config


from src.solver import Solver0D
from src.flow import Inflow0D
from src.misc import create_tool_parser, get_solver_path
from src.file_io import read_json
from src.stenosis import convert_all


@ray.remote
def proc_func(pbatch, changed_vessels, original_data, targets):
    results = []
    total = len(pbatch)
    
    # compute 
    #! not agnostic of where INFLOW is.
    tc = original_data['boundary_conditions'][0]['bc_values']['t'][-1] - original_data['boundary_conditions'][0]['bc_values']['t'][0]
    
    for pidx, p in enumerate(pbatch):
        # create original solver
        psolver = Solver0D()
        psolver.from_dict(copy.deepcopy(original_data))
        # fill psolver with correct modifications
        for idx, branch in enumerate(changed_vessels):
            for vessid in branch:
                vess = psolver.get_vessel(vessid)
                rad_rat = p[idx]
                convert_all(vess, rad_rat)
            
        print(f'Running sample {pidx + 1} / {total}')
        df = run_from_config(psolver.solver_data)
        
        
        # retrieve appropriate results
        

        
        internal_results = []
        
        for target_name in list(targets):
            temp_df = df[df['name'] == target_name]
            temp_df = temp_df[temp_df['time'] >= temp_df['time'].max() - tc]
            temp_df['time'] -= temp_df['time'].min()
            internal_results.append(np.trapz(temp_df['pressure_out'], temp_df['time']))
            internal_results.append(np.trapz(temp_df['flow_out'], temp_df['time']))
            
        results.append(internal_results)
        del psolver
    return pbatch, np.array(results)




    
def parametrize_stenosis(occlusion):
    ''' takes the artificial occlusion and computes the range of diameter changes for that particular vessel'''
    cur_area = 1 - occlusion
    area_increase = 1/cur_area
    radial_increase = np.sqrt(area_increase)
    # radial change goes from 1 (no repair) to (1/(1-occlusion))^(1/2) (complete repair)
    return 1, radial_increase


def data_gen(dims, occlusions, num_samples_log2 = 3):
    
    # generate sampler
    bits = 30 # means 2^30 max number of points... should be sufficient. Can go up to 64
    if num_samples_log2 > bits:
        raise ValueError(f'num_samples_log2 > {bits}, resulting in duplicates')
    sobol_sampler = Sobol(d= dims, scramble = True, bits = bits)
    
    # generate bounds
    lbound = []
    ubound = []
    for occ in occlusions:
        bounds = parametrize_stenosis(occ)
        lbound = bounds[0] - .05
        ubound = bounds[1]
    
    # create data
    # Make this a generator that does a max amount at a time
    return scale(sobol_sampler.random_base2(num_samples_log2), l_bounds=lbound, u_bounds=ubound)
        


        

    
    
def main(args):
    
    for solver_dir in tqdm.tqdm(args.solver_dirs, desc = 'Model'):
        
        join = lambda file: os.path.join(solver_dir, file)
        
        # get the solver/stenosis data
        solver_file = get_solver_path(solver_dir)
        solver_data = read_json(solver_file)
        stenosis_points = read_json(join("stenosis_vessels.dat"))
        
        # retrieve targets (last one of a branch)
        #! Subject to change
        targets = ['V' + str(x[-1]) for x in stenosis_points['all_changed_vessels']]


        #ray.init(num_cpus=args.num_proc)
        
        # get data
        data = data_gen(num_samples_log2=args.num_samples, dims = len(stenosis_points['all_changed_vessels']), occlusions = stenosis_points['occlusions'] )
        
        # divide data
        chunk_size = int(np.ceil(len(data) / args.num_proc))
        print(chunk_size)
        start = 0
        results = []
        
        for i in range(args.num_proc):
            # split data
            chunked_data = data[start:start + chunk_size]
            results.append(proc_func.remote(chunked_data, stenosis_points['all_changed_vessels'], solver_data, targets))
            start+=chunk_size
            
        results = ray.get(results)
        
        i = []
        o = []
        for input, output in results:
            i.append(input)
            o.append(output)

        if not os.path.exists(join('training_data')):
            os.mkdir(join('training_data'))
        
        np.save(join('training_data/input.npy'), np.array(np.concatenate(i, axis = 0)))
        np.save(join('training_data/output.npy'), np.array(np.concatenate(o, axis = 0)))
        
        
        ray.shutdown()
  

if __name__ == '__main__':
    
    parser = create_tool_parser(desc = 'Use Multi-threading to create data samples in numpy array for machine learning ')
    
    parser.add_argument('-n', dest = 'num_samples', default = 10, type = int, help = '2^num_samples will be generated. Default: 2^10 = 1024.')
    parser.add_argument('-p', dest = 'num_proc', default = 4, type = int, help = 'number of processes to use')
    
    args = parser.parse_args()
    main(args)