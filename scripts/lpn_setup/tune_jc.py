# File: tune_bc.py
# File Created: Monday, 31st October 2022 8:46:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 27th February 2023 11:48:38 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Tunes Boundary Conditions for a 0D model using a simplified tuning model.
#! Sensitivity Tests to be implemented & clean up code



from svinterface.core.zerod.solver import SolverResults, Solver0Dcpp
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.lpn import FastLPN, LPN
from svinterface.core.bc import Inflow, RCR
from svinterface.manager import Manager
from svinterface.utils.misc import m2d, d2m
from svinterface.utils.io import write_json

from pathlib import Path
import numpy as np
import argparse
import shutil
from copy import deepcopy
from scipy import optimize
import matplotlib.pyplot as plt

def opt_func(x, main_lpn: LPN, caps, targets):
    
    
    for junc in main_lpn.junctions:
        if junc['junction_type'] != 'internal_junction':
            
            S0 = junc['areas'][0]
            s_outlets = junc['areas'][1:]

            
            density = main_lpn.simulation_params['density']

            # update
            for idx, S1 in enumerate(s_outlets):
                jc = 1.52 * density * ( ((x[0] * S0)/(x[1] * S1) - 1) **2) / (2 * S0**2)
                main_lpn.change_junction_outlet(junc['junction_name'], which = idx, S = jc)
    
    s = Solver0Dcpp(main_lpn, last_cycle_only=True, mean_only=True)
    results = s.run_sim()
    values = [results.result_df.iloc[caps[0]]['pressure_in']]
    values += list(results.result_df.iloc[caps[1:]]['pressure_out'])
    values = np.array(values)
    loss = (((values - targets) / targets) ** 2).sum()
    print(x, loss)
    return loss


def tune(TM , main_lpn: LPN, threed_c: Centerlines):
    
    # get relevant positions
    tree = main_lpn.get_tree()
    # collect
    junction_vessels = []
    junction_gids = []
    junction_id = []
    for junc_node in main_lpn.tree_bfs_iterator(tree, allow='junction'):
        junction_vessels.append((junc_node.vessel_info[0]['inlet_vessels'], junc_node.vessel_info[0]['outlet_vessels']))
        junction_gids.append( junc_node.vessel_info[0]['gid']) # out gids
        junction_id.append(junc_node.id)
    
    
    # extract target resistances
    target_pressures = m2d(threed_c.get_pointdata_array("avg_pressure"))
    # get caps
    caps = []
    gids = []
    tree = main_lpn.get_tree()
    for node in main_lpn.tree_bfs_iterator(tree, allow = 'branch'):
        if 'boundary_conditions' in node.vessel_info[-1]:
            caps.append(node.ids[-1])
            if 'inlet' in node.vessel_info[-1]['boundary_conditions']:
                gids.append(node.vessel_info[-1]['gid'][0])
            else:
                gids.append(node.vessel_info[-1]['gid'][1])
    caps = np.array(caps)
    target_pressures = target_pressures[gids]
    #print(caps, target_pressures)
    bounds = optimize.Bounds([0,0], [5,5], keep_feasible=True)
    x0 = np.array([1,1])
    rez = optimize.minimize(opt_func,x0=x0, args =  (main_lpn, caps, target_pressures), bounds = bounds, method="Nelder-Mead",options={'disp':3} )
    
    print("X", rez.x, rez.fun)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Tunes a model if necessary')
    
    parser.add_argument('-i', dest = 'config', help = 'Config.yaml file')
    # dev params
    parser.add_argument('--f', dest = 'force', action = 'store_true', default = False, help = 'Whether to run an optimization even if tuning is set to false.')
    parser.add_argument('--s', dest = 'sensitivity_test', action = 'store_true', default = False, help = 'flag to run sensitivity tests or not')
    
    args = parser.parse_args()
    
    # manager
    TM = Manager(args.config)
    
    if args.sensitivity_test:
        raise NotImplementedError("Sensitivity Tests not implemented. Please remove the -s flag.")

        
    print("Tuning Model " + TM['metadata']['model_name'] + "...")
    
    # load the main LPN
    main_lpn = LPN.from_file(str(TM['workspace']['lpn']))
    # add arguments to params
    
    threed_file = TM['workspace']['3D']
    threed_c = Centerlines.load_centerlines(threed_file)
    # run optimizer
    tune(TM, main_lpn, threed_c)
