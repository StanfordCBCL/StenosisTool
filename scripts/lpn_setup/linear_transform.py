# File: linear_transform.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Sunday, 26th February 2023 4:35:27 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description:  Perform a linear transform on the junctions






import argparse

from svinterface.core.zerod.lpn import LPN
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
from svinterface.core.threed.extractor import Extract1D
import numpy as np
from concurrent.futures import ProcessPoolExecutor


def linear_transform(zerod_lpn: LPN, threed_c: Centerlines):

    # get relevant positions
    tree = zerod_lpn.get_tree()
    # collect
    junction_outlet_vessels = []
    junction_gids = []
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow='junction'):
        junction_outlet_vessels += junc_node.vessel_info[0]['outlet_vessels']
        junction_gids += junc_node.vessel_info[0]['gid'][1] # out gids
    
    assert len(junction_gids) == len(junction_outlet_vessels), "Number of junction ids on 3D data will not match the number of outlet vessels in the 0D"
        
    # extract target pressures.
    target_pressures = threed_c.get_pointdata_array("avg_pressure")[junction_gids]
    
    # compute initial case
    tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=True, debug = False)
    init_sim = tmp.run_sim()
    init_sim.convert_to_mmHg()
    pressures_init = init_sim.result_df.iloc[junction_outlet_vessels]['pressure_in'].to_numpy()
    
    
    pressures = []
    jcs = []
    # iterate through each junction outlet
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
        for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
            print(f"Changing junction {junc_node.id} vessel {idx}.")
            # compute a JC
            S0 = junc_node.vessel_info[0]['areas'][0]
            S1 = junc_node.vessel_info[0]['areas'][idx + 1]
            sc = 1.52 * zerod_lpn.simulation_params['density'] * ( (S0/S1 - 1) **2) / (2 * S0**2)
            jcs.append(sc)
            # change to JC
            zerod_lpn.change_junction_outlet(junction_id_or_name = junc_node.id, which = idx, S = sc)
            # Solver
            tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=True, debug = False)
            tmp_sim = tmp.run_sim()
            
            # get results
            tmp_sim.convert_to_mmHg()
            pressures_cur = tmp_sim.result_df.iloc[junction_outlet_vessels]['pressure_in'].to_numpy()
            pressures.append(pressures_cur - pressures_init) # append diff

            # undo change
            zerod_lpn.change_junction_outlet(junction_id_or_name = junc_node.id, which = idx, S = 0)
            
    # convert to numpy
    pressures = np.array(pressures).T
    
    # solve for a
    press_inv = np.linalg.inv(pressures) 
    
    # get target pressure differences
    target_pressures_diff = target_pressures - pressures_init
    
    aT = press_inv @ target_pressures_diff
    
    print(pressures)
    print(press_inv)
    print(target_pressures)
    print(aT)
    
    # iterate through each junction outlet and fill it with appropriate junction values
    counter = 0
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
        for idx in range(len(junc_node.vessel_info[0]['outlet_vessels'])):
            zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, S = aT[counter] * jcs[counter])
            counter += 1
            
    # save the lpn.
    zerod_lpn.update()
    
    

    

    
    
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Perform a linear optimization on the branches")
    parser.add_argument("-3d", dest = 'threed', help = '3d centerlines file')
    parser.add_argument("-0d", dest = 'zerod', help = "0d LPN")
    
    args = parser.parse_args()
    
    
    # get LPN and convert to BVJ
    zerod_lpn = LPN.from_file(args.zerod)
    zerod_lpn.to_cpp(normal = False) # resets to all 0's
    
    # load centerlines
    threed_c = Centerlines.load_centerlines(args.threed)
    
    linear_transform(zerod_lpn,threed_c)