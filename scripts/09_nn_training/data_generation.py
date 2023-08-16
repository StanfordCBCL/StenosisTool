# File: sobol_sampling_healthy.py
# File Created: Friday, 19th August 2022 4:22:32 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 15th August 2023 10:09:25 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Use Sobol sampling to parameterize from 0-1 each post-stent simulation. Save diastolic, mean, systolic pressures and flows.

import tqdm
import numpy as np
from scipy.stats.qmc import Sobol
from pathlib import Path
import json
import argparse
import math
from concurrent.futures import ProcessPoolExecutor


from svinterface.core.zerod.solver import Solver0Dcpp
from svinterface.core.zerod.lpn import LPN, FastLPN
from svinterface.manager.baseManager import Manager


def remote_run_sim(param, base_lpn: FastLPN, lpn_mapping: tuple, id: int):
    
    all_vess, all_vess_dr, all_juncs, all_juncs_dr = lpn_mapping
    all_targets = []
    # take the parameterization and apply it to lpn
    
    for pidx, p in enumerate(param):
        
        lpn = base_lpn.copy()
        print(f"Running simulation {pidx + 1}/{len(param)} on process {id}.", flush = True)
        # update lpns
        for idx, coef in enumerate(p):
            for vidx, max_dr in list(zip(all_vess[idx], all_vess_dr[idx])):
                # add the modified coefficient
                lpn.change_vessel(vidx, R = max_dr * coef, mode='add')

            for jidx, max_drs in list(zip(all_juncs[idx], all_juncs_dr[idx])):
                for outlet_idx, max_dr in enumerate(max_drs):
                    lpn.change_junction_outlet(int(jidx[1:]), which=outlet_idx, R=max_dr * coef, mode='add')
        
        
        # run simulation
        solver = Solver0Dcpp(lpn)
        results = solver.run_sim()
        
        # add MPA in
        vname = lpn.lpn_data[lpn.VESS][0]['vessel_name']
        targets = [*results.get_summ_val(vname, 'flow_in'),
                *results.get_summ_val(vname, 'pressure_in')]
        # for every vess, get the out
        for vess in lpn.lpn_data[base_lpn.VESS]:
            vname = vess['vessel_name']
            
            targets += [*results.get_summ_val(vname, 'flow_out')]
            targets += [*results.get_summ_val(vname, 'pressure_out')]
            
        all_targets.append(targets)

        del lpn
        
    return np.array(all_targets)

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
        # load the regions used
        rrpath = Path(M['parameterization']['corrections'][sim]['relevant_regions'])
        with rrpath.open() as rrfp:
            rr = json.load(rrfp)
        # load lpn
        lpn = LPN.from_file(M['parameterization']['corrections'][sim]['lpn'])
        all_vess.append(rr['Vessels'])
        all_vess_dr.append([lpn.get_vessel(vidx)['zero_d_element_values']['R_poiseuille'] - base_lpn.get_vessel(vidx)['zero_d_element_values']['R_poiseuille'] for vidx in rr['Vessels']])
        
        all_juncs.append(rr['Junctions'])
        all_juncs_dr.append([[lpn.get_junction(jidx)['junction_values']['R_poiseuille'][outlet_idx] - base_lpn.get_junction(jidx)['junction_values']['R_poiseuille'][outlet_idx] for outlet_idx in range(len(lpn.get_junction(jidx)['outlet_vessels']))] for jidx in rr['Junctions']])
            
    return base_lpn, all_vess, all_vess_dr, all_juncs, all_juncs_dr
            


def generate_data(M: Manager, data_dir: Path, samples: list, nprocs: int):
    
    train_dir = data_dir / 'train_data'
    val_dir = data_dir / 'val_data'
    test_dir = data_dir / 'test_data'
    
    total_sims = len(get_sim_names(M))
    base_lpn, all_vess, all_vess_dr, all_juncs, all_juncs_dr = parameterize(M)
    
    for idx, (mode_dir, num_samples) in enumerate(zip([train_dir, val_dir, test_dir], samples)):
        
        # make dir
        mode_dir.mkdir(exist_ok=True)
        
        # sobol sample data
        parameterization = sobol_data_gen(size=total_sims,
                              num_samples=num_samples,
                              seed=42 + idx)

        # pass to each process, each process handles incr (64) simulations before creating a fresh process

        # split data
        chunked_data = [[] for i in range(nprocs)]
        for idx, d in enumerate(parameterization):
            chunked_data[idx % nprocs].append(d)

        # pass to each process
        with ProcessPoolExecutor(max_workers=nprocs) as executor:
            futures = executor.map(remote_run_sim, chunked_data, [base_lpn.get_fast_lpn() for i in range(nprocs)], [(all_vess, all_vess_dr, all_juncs, all_juncs_dr) for i in range(nprocs)], range(nprocs))

        y = []
        # get y's
        for f in futures:
            y.append(f.result())
            
        y = np.concatenate(y)
        np.save(mode_dir / 'input.npy', parameterization)
        np.save(mode_dir / 'output.npy', y)
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Use Multi-threading to create data samples in numpy array for machine learning ')
    
    parser.add_argument('-i', dest = 'config', required = True, help = 'yaml file for confit')
    parser.add_argument('-ntrain', dest = 'num_train_samples', default = 8192, type = int, help = 'num_train_samples will be generated for training data. Use a power of 2 to guarentee balance properties. Default: 8192 = 2^13.')
    parser.add_argument('-nval', dest = 'num_val_samples', default = 1024, type = int, help = 'num_val_samples will be generated for validation data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-ntest', dest = 'num_test_samples', default = 1024, type = int, help = 'num_test_samples will be generated for testing data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-nprocs', type=int, required=True, help = 'number of processes to use')
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    nn_dir = Path(M['workspace']['nn_dir'])
    
    # main data dir
    data_dir = nn_dir / 'model_data'
    data_dir.mkdir(exist_ok=True)
    M.register('model_data', str(data_dir), depth=['NN_DIR'])
    
    generate_data(M, data_dir, [args.num_train_samples, args.num_val_samples, args.num_test_samples], args.nprocs)
    
    