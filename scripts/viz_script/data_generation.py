# File: sobol_sampling_healthy.py
# File Created: Friday, 19th August 2022 4:22:32 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 25th September 2023 3:04:07 pm
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


def change_sim(param, base_lpn: FastLPN, lpn_mapping: tuple):
    all_vess, all_vess_dr, all_juncs, all_juncs_dr = lpn_mapping
    # take the parameterization and apply it to lpn

    for idx, coef in enumerate(param):
        for vidx, max_dr in list(zip(all_vess[idx], all_vess_dr[idx])):
            # add the modified coefficient
            base_lpn.change_vessel(vidx, R = max_dr * coef, mode='add')

        for jidx, max_drs in list(zip(all_juncs[idx], all_juncs_dr[idx])):
            for outlet_idx, max_dr in enumerate(max_drs):
                base_lpn.change_junction_outlet(int(jidx[1:]), which=outlet_idx, R=max_dr * coef, mode='add')
                
    return base_lpn

def get_sim_names(M: Manager):
    ''' gets names of simulations'''
    with open("data/diseased/AS1_SU0308_stent/results/AS1_SU0308_nonlinear/NN_DIR_old/parametrization_sims.txt", 'r') as sfile:
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
        print(all_juncs)
        all_juncs_dr.append([[lpn.get_junction(jidx)['junction_values']['R_poiseuille'][outlet_idx] - base_lpn.get_junction(jidx)['junction_values']['R_poiseuille'][outlet_idx] for outlet_idx in range(len(lpn.get_junction(jidx)['outlet_vessels']))] for jidx in rr['Junctions']])
        print(all_juncs_dr)
    return base_lpn, all_vess, all_vess_dr, all_juncs, all_juncs_dr
            


def generate_data(M: Manager):
    """ Generate Data Proc """
    
    base_lpn, all_vess, all_vess_dr, all_juncs, all_juncs_dr = parameterize(M)
    # sobol sample data
    parameterization = np.array([1.0325e-01, 2.1508e-01, 5.3708e-01])
    # pass to each process, each executor handles incr (128) simulations before creating a fresh set of subprocesses

    base_lpn = change_sim(parameterization, base_lpn, (all_vess, all_vess_dr, all_juncs, all_juncs_dr))
    
    base_lpn.write_lpn_file("test.in")
    
    
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Use Multi-threading to create data samples in numpy array for machine learning ')
    parser.add_argument('-i', dest = 'config', required = True, help = 'yaml file for confit')
    parser.add_argument('-ntrain', dest = 'num_train_samples', default = 8192, type = int, help = 'num_train_samples will be generated for training data. Use a power of 2 to guarentee balance properties. Default: 8192 = 2^13.')
    parser.add_argument('-nval', dest = 'num_val_samples', default = 1024, type = int, help = 'num_val_samples will be generated for validation data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    parser.add_argument('-ntest', dest = 'num_test_samples', default = 1024, type = int, help = 'num_test_samples will be generated for testing data. Use a power of 2 to guarentee balance properties. Default: 1024 = 2^10.')
    args = parser.parse_args()
    
    M = Manager(args.config)


    # generate data
    generate_data(M)
    
    