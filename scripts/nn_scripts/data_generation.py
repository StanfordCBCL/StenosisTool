# File: sobol_sampling_healthy.py
# File Created: Friday, 19th August 2022 4:22:32 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Saturday, 3rd December 2022 2:34:10 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Use Sobol sampling to retrieve a distribution of potential diameter changes for a particular healthy model. Takes in an artificial stenosis directory.

import tqdm
import ray
import os
import pandas as pd
import numpy as np
from scipy.stats.qmc import Sobol, scale
import copy
from pathlib import Path


from sgt.core.solver import Solver0Dcpp, SolverResults
from sgt.core.lpn import LPN
from sgt.utils.parser import ToolParser
from sgt.core.manager import DataGenerationManager
from sgt.utils.io import read_json

@ray.remote
def proc_func(pbatch, changed_vessels, original_data: LPN, targets):
    results = []
    total = len(pbatch)
    
    # compute
    for pidx, p in enumerate(pbatch):
        # create copy of original solver
        pLPN = original_data.deep_copy()
        # fill psolver with correct modifications
        for idx, branch in enumerate(changed_vessels):
            for vidx in branch:
                pLPN.repair_vessel(vidx, p[idx]**2 ) # repair by rad^2
            
        print(f'Running sample {pidx + 1} / {total}', flush = True)
        # solve
        solver = Solver0Dcpp(lpn=pLPN,
                             last_cycle_only=True,
                             mean_only=True,
                             debug=False)
        solver_results=solver.run_sim()
        

        # retrieve results
        internal_results = []
        
        #! Modify this depending on what is desired (Currently only mean location after)
        for (in_target, out_target) in list(targets):
            temp_df = solver_results.vessel_df(out_target)
            internal_results.append(temp_df['pressure_out'][0])
            internal_results.append(temp_df['flow_out'][0])
            
            
            
        results.append(internal_results)
        del pLPN
    return pbatch, np.array(results)



def data_gen(dims, max_repairs, num_samples, seed, lower_offset = 0.05):
    
    # generate sampler
    bits = 30 # means 2^30 max number of points... should be sufficient. Can go up to 64
    if num_samples > 2**bits:
        raise ValueError(f'num_samples > {2**bits}, resulting in duplicates')
    sobol_sampler = Sobol(d= dims, scramble = True, bits = bits, seed = seed)
    
    lbound = [1 - lower_offset for i in range(len(max_repairs))]
    ubound = max_repairs
    
    # create data
    return scale(sobol_sampler.random(num_samples), l_bounds=lbound, u_bounds=ubound)

def get_targets(vessels):
    ''' gets various targets '''
    in_targets = []
    out_targets = []
    for v in vessels:
        in_targets.append(data_lpn.get_vessel(v[0])['vessel_name'])
        out_targets.append(data_lpn.get_vessel(v[-1])['vessel_name'])
        
    return zip(in_targets, out_targets)

def parse_parametrization(parametrization):
    ''' parses the parametrization file '''
    max_repair = []
    vessels = []
    for i in range( len(parametrization) ):
        p = parametrization[str(i)]
        # get the targets for each vessel to save
        
        # get vess id
        vessels.append(p[0])
        # get max repair
        max_repair.append(p[1])

    return  vessels,max_repair
    

if __name__ == '__main__':
    
    parser = ToolParser(desc = 'Use Multi-threading to create data samples in numpy array for machine learning ')
    
    parser.parser.add_argument('-ntrain', dest = 'num_train_samples', default = 1024, type = int, help = 'num_train_samples will be generated for training data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.parser.add_argument('-nval', dest = 'num_val_samples', default = 1024, type = int, help = 'num_val_samples will be generated for validation data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.parser.add_argument('-ntest', dest = 'num_test_samples', default = 1024, type = int, help = 'num_test_samples will be generated for testing data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.parser.add_argument('-p', dest = 'num_proc', default = 4, type = int, help = 'number of processes to use')
    
    args = parser.parse_args()
    
    M = DataGenerationManager(args.config)
    
    # LPN
    data_lpn = LPN.from_file(str(M.data_gen_lpn))
    
    # parametrization
    parametrization = read_json(str(M.stenosis_parametrization))
    
    vessels, max_repairs = parse_parametrization(parametrization)
    
    targets = get_targets(vessels)
    
    # start ray
    ray.init(num_cpus=args.num_proc)
        
    for idx, (mode, num_samples) in enumerate([('train', args.num_train_samples), ('val',args.num_val_samples), ('test', args.num_test_samples)]):
            
            mode_dir = M.data_dir / mode
            mode_dir.mkdir(exists_ok = True)
            
            # change offset 
            if mode == 'train':
                lower_offset = 0.05
            else:
                lower_offset = 0
                
            # get data (with different seeds)
            data = data_gen(num_samples=num_samples, 
                            dims = len(parametrization), 
                            max_repairs = max_repairs,
                            seed = 42 + idx,
                            lower_offset=lower_offset)
            
            # divide data
            chunk_size = int(np.ceil(len(data) / args.num_proc))
            start = 0
            results = []
            
            for i in range(args.num_proc):
                # split data
                chunked_data = data[start:start + chunk_size]
                results.append(proc_func.remote(chunked_data, vessels, data_lpn, targets))
                start+=chunk_size
            results = ray.get(results)
            
            i = []
            o = []
            for input, output in results:
                i.append(input)
                o.append(output)


            
            np.save(mode_dir / 'input.npy', np.array(np.concatenate(i, axis = 0)))
            np.save(mode_dir / 'output.npy', np.array(np.concatenate(o, axis = 0)))
            
    ray.shutdown()