# File: linear_transform.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 27th February 2023 2:31:21 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description:  Perform a linear transform on the junctions






import argparse

from svinterface.core.zerod.lpn import LPN
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
from svinterface.manager.baseManager import Manager
from svinterface.utils.misc import d2m, m2d
from svinterface.core.threed.extractor import Extract1D
import numpy as np
from concurrent.futures import ProcessPoolExecutor


def linear_transform(zerod_lpn: LPN, threed_c: Centerlines):

    # get relevant positions
    tree = zerod_lpn.get_tree()
    # collect
    junction_vessels = []
    junction_gids = []
    junction_id = []
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow='junction'):
        junction_vessels.append((junc_node.vessel_info[0]['inlet_vessels'], junc_node.vessel_info[0]['outlet_vessels']))
        junction_gids.append( junc_node.vessel_info[0]['gid']) # out gids
        junction_id.append(junc_node.id)

    
    # Solver
    tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=True, debug = False)
    tmp_sim = tmp.run_sim()
    mPAP = tmp_sim.result_df.iloc[0]['pressure_in']
    
    
    # extract target resistances
    target_pressures = m2d(threed_c.get_pointdata_array("avg_pressure"))
    pressures_0d = tmp_sim.result_df['pressure_in']
    pressures_0d_out = tmp_sim.result_df['pressure_out']
    r_diff = (target_pressures[0] - mPAP)/ tmp_sim.result_df.iloc[0]['flow_in']
    
    for idx, (in_gid, out_gid) in enumerate(junction_gids):
        print(target_pressures[out_gid], target_pressures[in_gid])
        dP =  -(target_pressures[out_gid] - target_pressures[in_gid])
        Q = tmp_sim.result_df.iloc[junction_vessels[idx][1]]['flow_in'].to_numpy()
        
        # compute R's
        R = dP / Q
        j = zerod_lpn.get_junction(junction_id[idx])
        for o_idx,_ in enumerate(out_gid):
            zerod_lpn.change_junction_outlet(junction_id_or_name = junction_id[idx], which = o_idx, R = R[o_idx])
    
    zerod_lpn.update()
    
    exit(1)
    
    pressures = []
    jcs = []
    # iterate through each junction outlet
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
        for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
            print(f"Changing junction {junc_node.id} vessel {idx}.")
            # compute a JC
            r = 1
            jcs.append(r)
            # change to JC
            zerod_lpn.change_junction_outlet(junction_id_or_name = junc_node.id, which = idx, R = r)
            # Solver
            tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=True, debug = False)
            tmp_sim = tmp.run_sim()
            
            # get results
            tmp_sim.convert_to_mmHg()
            pressures_cur = tmp_sim.result_df.iloc[junction_outlet_vessels]['pressure_in'].to_numpy()
            pressures.append(pressures_cur - pressures_init) # append diff

            # undo change
            zerod_lpn.change_junction_outlet(junction_id_or_name = junc_node.id, which = idx, R = 0)
            
    # convert to numpy
    pressures = np.array(pressures).T
    
    # solve for a
    press_inv = np.linalg.inv(pressures) 
    
    # get target pressure differences
    target_pressures_diff = target_pressures - pressures_init
    
    aT = press_inv @ target_pressures_diff
    
    print(target_pressures_diff)
    print(aT)
    
    # iterate through each junction outlet and fill it with appropriate junction values
    counter = 0
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
        for idx in range(len(junc_node.vessel_info[0]['outlet_vessels'])):
            zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, R = aT[counter] * jcs[counter])
            counter += 1
            
    # save the lpn.
    zerod_lpn.update()
    
    

    

    
    
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Perform a linear optimization on the branches")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    
    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    zerod_file = M['workspace']['lpn']
    threed_file = M['workspace']['3D']
    
    # get LPN and convert to BVJ
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.to_cpp(normal = False) # resets to all 0's
    zerod_lpn.update()
    
    # load centerlines
    threed_c = Centerlines.load_centerlines(threed_file)
    linear_transform(zerod_lpn,threed_c)