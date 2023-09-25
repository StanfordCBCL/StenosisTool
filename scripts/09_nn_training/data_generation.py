# File: sobol_sampling_healthy.py
# File Created: Friday, 19th August 2022 4:22:32 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 25th September 2023 9:07:04 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Use Sobol sampling to parameterize from 0-1 each post-stent simulation. Save diastolic, mean, systolic pressures and flows.

import numpy as np
from scipy.stats.qmc import Sobol
from pathlib import Path
import json
import argparse
import time
from concurrent.futures import ProcessPoolExecutor

from svinterface.core.zerod.solver import Solver0Dcpp
from svinterface.core.zerod.lpn import LPN, FastLPN
from svinterface.manager.baseManager import Manager


def remote_run_sim(param, base_lpn: FastLPN, lpn_mapping: tuple):
    
    all_vess, all_vess_dr, all_juncs, all_juncs_dr = lpn_mapping
    # take the parameterization and apply it to lpn
    
    for idx, coef in enumerate(param):
        for vidx, max_dr in list(zip(all_vess[idx], all_vess_dr[idx])):
            # add the modified coefficient
            base_lpn.change_vessel(vidx, R = max_dr * coef, mode='add')
    
        for jidx, max_drs in list(zip(all_juncs[idx], all_juncs_dr[idx])):
            for outlet_idx, max_dr in enumerate(max_drs):
                base_lpn.change_junction_outlet(int(jidx[1:]), which=all_juncs[idx][jidx][outlet_idx], R=max_dr * coef, mode='add')
    
    
    # run simulation
    solver = Solver0Dcpp(base_lpn)
    results = solver.run_sim()
    
    # add MPA in
    vname = base_lpn.lpn_data[base_lpn.VESS][0]['vessel_name']
    targets = [*results.get_summ_val(vname, 'flow_in'),
            *results.get_summ_val(vname, 'pressure_in')]
    # for every vess, get the out
    for vess in base_lpn.lpn_data[base_lpn.VESS]:
        vname = vess['vessel_name']
        
        targets += [*results.get_summ_val(vname, 'flow_out')]
        targets += [*results.get_summ_val(vname, 'pressure_out')]
    
    del results, solver

    return np.float32(np.array(targets))

def sobol_data_gen(size, num_samples, seed):
    
    # generate sampler
    bits = 30 # means 2^30 max number of points... should be sufficient. Can go up to 64
    if num_samples > 2**bits:
        raise ValueError(f'num_samples > {2**bits}, resulting in duplicates')
    sobol_sampler = Sobol(d=size, scramble=True, bits=bits, seed=seed)
    
    # sample
    return sobol_sampler.random(num_samples)

def get_sim_names(M: Manager):
    ''' gets names of simulations'''
    with open(M['NN_DIR']['sims'], 'r') as sfile:
        sims = sfile.readlines()
    return [s.rstrip() for s in sims]

def parameterize(M: Manager):
    sim_names = get_sim_names(M)
    
    # load base_lpn
    base_lpn = LPN.from_file(M['parameterization']['base_lpn'])

    all_vess = []
    all_vess_dr = []
    all_juncs = []
    all_juncs_dr = []
    # pull out relevant regions and get vessel and junction resistance ranges
    for sim in sim_names:
        print(sim_names)
        # load the regions used
        rrpath = Path(M['parameterization']['corrections'][sim]['relevant_regions'])
        with rrpath.open() as rrfp:
            rr = json.load(rrfp)
            
        print(rr)
        # load lpn
        lpn = LPN.from_file(M['parameterization']['corrections'][sim]['lpn'])
        all_vess.append(rr['Vessels'])
        all_vess_dr.append([lpn.get_vessel(vidx)['zero_d_element_values']['R_poiseuille'] - base_lpn.get_vessel(vidx)['zero_d_element_values']['R_poiseuille'] for vidx in rr['Vessels']])
        
        all_juncs.append(rr['Junctions'])
        all_juncs_dr.append([[lpn.get_junction(jidx)['junction_values']['R_poiseuille'][outlet_idx] - base_lpn.get_junction(jidx)['junction_values']['R_poiseuille'][outlet_idx] for outlet_idx in rr['Junctions'][jidx]] for jidx in rr['Junctions']])
        print(all_juncs)
    return base_lpn, all_vess, all_vess_dr, all_juncs, all_juncs_dr
            


def generate_data(M: Manager, data_dir: Path, samples: list):
    """ Generate Data Proc """
    train_dir = data_dir / 'train_data'
    val_dir = data_dir / 'val_data'
    test_dir = data_dir / 'test_data'
    
    total_sims = len(get_sim_names(M))
    base_lpn, all_vess, all_vess_dr, all_juncs, all_juncs_dr = parameterize(M)
    
    for idx, (name, mode_dir, num_samples) in enumerate(zip(['train data', 'val data', 'test data'],[train_dir, val_dir, test_dir], samples)):
        
        # make dir
        mode_dir.mkdir(exist_ok=True)
        
        # sobol sample data
        parameterization = sobol_data_gen(size=total_sims,
                              num_samples=num_samples,
                              seed=42 + idx)

        # pass to each process, each executor handles incr (128) simulations before creating a fresh set of subprocesses

        y = []
        cur, incr = 0, 32
        counter = 0
        
        while cur < num_samples:
            
            # restart executor
            if cur % (4 * incr) == 0:
                executor = ProcessPoolExecutor()

            # submit jobs in batches of 32
            start = time.time()
            futures = []
            for p in parameterization[cur:cur + incr]:
                futures.append(executor.submit(remote_run_sim, p, base_lpn.get_fast_lpn(), (all_vess, all_vess_dr, all_juncs, all_juncs_dr)))
            
            # get futures
            for f in futures:
                y.append(f.result())
                counter += 1
                print(f"Retrieved results for simulation {counter}/{num_samples}.", flush = True)
                
            
            del futures
            cur += incr
            print(f"Time (sec) taken for last {incr} jobs: {time.time() - start}")
            
            # restart executor
            if cur % (4 * incr) == 0:
                executor.shutdown(wait=True)
            
        y = np.vstack(y)
        np.save(mode_dir / 'input.npy', parameterization)
        np.save(mode_dir / 'output.npy', y)
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Use Multi-threading to create data samples in numpy array for machine learning ')
    parser.add_argument('-i', dest = 'config', required = True, help = 'yaml file for confit')
    parser.add_argument('-ntrain', dest = 'num_train_samples', default = 8192, type = int, help = 'num_train_samples will be generated for training data. Use a power of 2 to guarentee balance properties. Default: 8192 = 2^13.')
    parser.add_argument('-nval', dest = 'num_val_samples', default = 1024, type = int, help = 'num_val_samples will be generated for validation data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-ntest', dest = 'num_test_samples', default = 1024, type = int, help = 'num_test_samples will be generated for testing data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    # main data dir
    nn_dir = Path(M['workspace']['nn_dir'])
    data_dir = nn_dir / 'model_data'
    data_dir.mkdir(exist_ok=True)
    M.register('model_data', str(data_dir), depth=['NN_DIR'])
    
    # generate data
    generate_data(M, data_dir, [args.num_train_samples, args.num_val_samples, args.num_test_samples])
    
    