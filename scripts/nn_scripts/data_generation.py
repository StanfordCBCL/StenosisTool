# File: sobol_sampling_healthy.py
# File Created: Friday, 19th August 2022 4:22:32 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 3rd April 2023 9:38:33 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Use Sobol sampling to retrieve a distribution of potential diameter changes for a particular healthy model. Takes in an artificial stenosis directory.

import tqdm
import concurrent.futures 
import os
import pandas as pd
import numpy as np
from scipy.stats.qmc import Sobol, scale
import copy
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor


from svinterface.core.zerod.solver import Solver0Dcpp
from svinterface.core.zerod.lpn import LPN, FastLPN
from svinterface.manager.baseManager import Manager
from svinterface.utils.io import read_json


def proc_func(data, changed_vessels, lpn: FastLPN, targets):
    results = []
    # compute
    for didx, d in enumerate(data):
        # create copy of original solver
        tmpLPN = lpn.copy()
        # fill psolver with correct modifications
        for idx, branch in enumerate(changed_vessels):
            for vidx in branch:
                tmpLPN.repair_vessel(vidx, d[idx]**2) # repair by rad^2
            
        print(f'Running sample {didx + 1} / {len(data)}', flush = True)
        # solve
        solver = Solver0Dcpp(lpn=tmpLPN,
                             last_cycle_only=True,
                             mean_only=True,
                             debug=False)
        solver_results=solver.run_sim()
        

        # retrieve results
        in_target, out_target = targets
        internal_results = []
        
        #! Modify this depending on what is desired (Currently only mean location after)
        for vess_list in in_target:
            results.append(solver_results.result_df.iloc[vess_list[0]]['pressure_in'])
            results.append(solver_results.result_df.iloc[vess_list[0]]['flow_in'])
            
        for vess_list in out_target:
            results.append(solver_results.result_df.iloc[vess_list[-1]]['pressure_out'])
            results.append(solver_results.result_df.iloc[vess_list[-1]]['flow_out'])
                 
        results.append(internal_results)
        del tmpLPN
    return np.array(results)



def data_gen( max_repairs, num_samples, seed, lower_offset = 0.05):
    
    # generate sampler
    bits = 30 # means 2^30 max number of points... should be sufficient. Can go up to 64
    if num_samples > 2**bits:
        raise ValueError(f'num_samples > {2**bits}, resulting in duplicates')
    sobol_sampler = Sobol(d= len(max_repairs), scramble = True, bits = bits, seed = seed)
    
    lbound = np.ones(len(max_repairs)) - lower_offset
    ubound = max_repairs
    
    # create data
    return scale(sobol_sampler.random(num_samples), l_bounds=lbound, u_bounds=ubound)
    
def main(args):
    
    M = Manager(args.config)
    
    
    if args.mode == 'AS':
        sims = 'as_simulations'
    elif args.mode == 'R':
        sims = 'r_simulations'
    else:
        raise ValueError("-mode must be AS or R")

    S = M[sims][args.sim]
    register_main = lambda key, val: M.register(key, val, depth = [sims, args.sim])
    
    # dir
    main_dir = Path(S['dir'])
    # lpn
    lpn = LPN.from_file(S['lpn'])
    
    ## read param and parse into two arrays
    stenp = read_json(Path(S['sten_p']))
    vessels, max_repairs = zip(*stenp)
    vessels, max_repairs = list(vessels), list(max_repairs)
    
    
    ## get targets
    in_targets = [[0]]
    out_targets = vessels
    
    
    ## create directories for data
    register_main('model_data', {})
    
    # main data dir
    data_dir = main_dir / 'model_data'
    data_dir.mkdir(exist_ok=True)
    register_data = lambda key, val: M.register(key, val, depth = [sims, args.sim, 'model_data'])
    register_data('dir', str(data_dir))    

    # subdirs
    train_dir = data_dir / 'train_data'
    val_dir = data_dir / 'val_data'
    test_dir = data_dir / 'test_data'
    
    for idx, (mode, mode_dir, num_samples ) in enumerate([('train',train_dir, args.num_train_samples), ('val', val_dir, args.num_val_samples), ('test', test_dir, args.num_test_samples)]):
        
        # make dir
        mode_dir.mkdir(exist_ok=True)
        
        # compute an offset for train
        if mode == 'train':
            lower_offset = 0.05
        else:
            lower_offset = 0
        
        # sobol sample data
        data = data_gen(num_samples=num_samples, 
                            max_repairs = max_repairs,
                            seed = 42 + idx,
                            lower_offset=lower_offset)
        
        # split data
        chunked_data = [[] for i in range(args.p)]
        for idx, d in enumerate(data):
            chunked_data[idx % args.p].append(d)

        # pass to each process
        with ProcessPoolExecutor(max_workers=args.p) as executor:
            futures = executor.map(proc_func, chunked_data, [vessels for i in range(args.p)], [lpn.get_fast_lpn() for i in range(args.p)], [(in_targets, out_targets) for i in range(args.p)] )
        
        # get x, y's
        x = np.concatenate(chunked_data)
        y = []
        for f in futures:
            y.append(f.result())
        y = np.concatenate(y)
        
        np.save(mode_dir / 'input.npy', x)
        np.save(mode_dir / 'output.npy', y)
        
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Use Multi-threading to create data samples in numpy array for machine learning ')
    
    parser.add_argument('-i', dest = 'config', required = True, help = 'yaml file for confit')
    parser.add_argument('-mode', required=True, help = 'Mode: AS, R')
    parser.add_argument("-sim", type = int, help = "Simulation number to use")
    parser.add_argument('-ntrain', dest = 'num_train_samples', default = 1024, type = int, help = 'num_train_samples will be generated for training data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-nval', dest = 'num_val_samples', default = 1024, type = int, help = 'num_val_samples will be generated for validation data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-ntest', dest = 'num_test_samples', default = 1024, type = int, help = 'num_test_samples will be generated for testing data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-p', default = 4, type = int, help = 'number of processes to use')
    
    args = parser.parse_args()
    
    main(args)
    