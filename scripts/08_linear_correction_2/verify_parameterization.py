# File: parametrization.py
# File Created: Friday, 4th August 2023 12:04:58 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 14th August 2023 5:16:07 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Verifies that the parameterization used to construct a neural network is completely valid


import argparse
import warnings
from pathlib import Path
import json

from svinterface.manager import Manager
from svinterface.core.zerod.lpn import LPN


def check_overlapping_regions(M: Manager, sims: list):
    
    ## first find the regions
    overlap_vess, overlap_juncs = set(), set()
    all_vess, all_juncs = set(), set()
    for sim in sims:
        # load the regions used
        rrpath = Path(M['parameterization']['corrections'][sim]['relevant_regions'])
        with rrpath.open() as rrfp:
            rel_regions = json.load(rrfp)
        overlap_vess = overlap_vess.union(all_vess.intersection(set(rel_regions['Vessels'])))
        all_vess = all_vess.union(set(rel_regions['Vessels']))
        overlap_juncs = overlap_juncs.union(all_juncs.intersection(set(rel_regions['Junctions'])))
        all_juncs = all_juncs.union(set(rel_regions['Junctions']))

    if len(overlap_vess) > 0:
        print(f"Warning: there are overlapping vessels at the following regions: {overlap_vess}")
    
    if len(overlap_juncs) > 0:
        print(f"Warning: there are overlapping junctions at the following regions: {overlap_juncs}")

def check_param(M, sims):
    """ checks minimum possible resistance achievable """
    base_lpn = LPN.from_file(M['parameterization']['base_lpn'])
    lpns = []
    for sim in sims:
        lpns.append(LPN.from_file(M['parameterization']['corrections'][sim]['lpn']))
    
    get_r = lambda vess: vess['zero_d_element_values']['R_poiseuille']
    get_rj = lambda junc, idx: junc['junction_values']['R_poiseuille'][idx]
    for vess_idx in range(base_lpn.num_vessels()):
        # check if all the negative/0 terms linearly combined can reach a negative resistance
        if (get_r(base_lpn.get_vessel(vess_idx)) + sum([min(0, get_r(lpn.get_vessel(vess_idx)) - get_r(base_lpn.get_vessel(vess_idx))) for lpn in lpns])) < 0 :
            raise ValueError(f"At vessel {vess_idx}, an invalid resistance may occur")
    
    for junc_idx in range(base_lpn.num_junctions()):
        # skip internal junctions
        if base_lpn.get_junction(junc_idx)['junction_type'] == 'internal_junction':
            continue
        
        # check if all the negative/0 terms linearly combined can reach a negative resistance
        for outlet_idx in range(len(base_lpn.get_junction(junc_idx)['outlet_vessels'])):
            if (get_rj(base_lpn.get_junction(junc_idx), outlet_idx) + sum([min(0, get_rj(lpn.get_junction(junc_idx), outlet_idx) - get_rj(base_lpn.get_junction(junc_idx), outlet_idx)) for lpn in lpns])) < 0 :
                raise ValueError(f"At junction J{junc_idx} outlet {outlet_idx}, an invalid resistance may occur")
            
        
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Parametrizes the diseased to each stented as LPN")
    
    parser.add_argument("-i", dest='config', help = 'yaml config file')
    parser.add_argument("-sims", nargs='+', help = 'name of simulations from parameterization to use in model')
    args = parser.parse_args()
    
    M=Manager(args.config)
    
    # add nndir cfg
    if 'NN_DIR' not in M.yaml or M['NN_DIR'] is None:
        M.register('NN_DIR', {})
    
    # add nn dir
    nn_dir = Path(M['workspace']['root']) / 'NN_DIR'
    nn_dir.mkdir(exist_ok=True)
    M.register('nn_dir', str(nn_dir), ['workspace'])
    
    # check used parameterizations
    check_overlapping_regions(M, args.sims)
    check_param(M, args.sims)
    
    # add the relevant sims used
    param_sims = nn_dir / 'parametrization_sims.txt'
    with param_sims.open("w") as pfile:
        pfile.write("\n".join(list(args.sims)))
    M.register('sims', str(param_sims), ['NN_DIR'])
        
        
    print("Parameterization Verified!!!")
    
    M.update()
    
    