# File: linear_transform.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 29th June 2023 1:22:46 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description:  Perform a linear transform on the junctions






import argparse

from svinterface.core.zerod.lpn import LPN, OriginalLPN
from svinterface.core.polydata import Centerlines
from svinterface.core.bc import RCR 
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
from svinterface.manager.baseManager import Manager
from svinterface.utils.misc import m2d
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, wait



def junc_sim(lpn: OriginalLPN, dvess: int, junc_id: int, which: int, poi_0d):
    ''' To run simulations using multi-processing '''
    # get downstream vessel.
    v = lpn.get_vessel(dvess)
    # compute a R as 10%
    dr = .1 * v['zero_d_element_values']['R_poiseuille']
    # change to JC
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = dr, mode = 'add')
    # Solver
    tmp = Solver0Dcpp(lpn, last_cycle_only=True, mean_only=True, debug = False)
    tmp_sim = tmp.run_sim()

    # get results
    tmp_sim.convert_to_mmHg()
    pressures_cur = pd.concat([tmp_sim.result_df.iloc[poi_0d[0]]['pressure_in'], tmp_sim.result_df.iloc[poi_0d[1]]['pressure_out'], tmp_sim.result_df.iloc[poi_0d[2]]['pressure_in']]).to_numpy()
    # undo change
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = -dr, mode = 'add')

    return dr, pressures_cur   

def vess_sim(lpn: OriginalLPN, vess: int, poi_0d):
    ''' To run simulations using multi-processing '''
    # get downstream vessel.
    v = lpn.get_vessel(vess)
    # compute a R as 10%
    dr = .1 * v['zero_d_element_values']['R_poiseuille']
    # change to JC
    lpn.change_vessel(vessel_id = vess, R = dr, mode = 'add')
    # Solver
    tmp = Solver0Dcpp(lpn, last_cycle_only=True, mean_only=True, debug = False)
    tmp_sim = tmp.run_sim()

    # get results
    tmp_sim.convert_to_mmHg()
    pressures_cur = pd.concat([tmp_sim.result_df.iloc[poi_0d[0]]['pressure_in'], tmp_sim.result_df.iloc[poi_0d[1]]['pressure_out'], tmp_sim.result_df.iloc[poi_0d[2]]['pressure_in']]).to_numpy()
    # undo change
    lpn.change_vessel(vessel_id = vess, R = -dr, mode = 'add')

    return dr, pressures_cur       
        
def linear_transform(zerod_lpn: LPN, threed_c: Centerlines, M: Manager):

    # get relevant positions
    tree = zerod_lpn.get_tree()
    # collect junction values
    junction_outlet_vessels = []
    junction_gids = []
    junction_nodes = []
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow='junction'):
        junction_outlet_vessels += junc_node.vessel_info[0]['outlet_vessels']
        junction_gids += junc_node.vessel_info[0]['gid'][1] # out gids
        junction_nodes.append(junc_node)
        
    # collect branch values
    segment_vess_ids = []
    segment_gids = []
    branch_nodes = []
    for branch_node in zerod_lpn.tree_bfs_iterator(tree, allow='branch'):
        
        for idx, vess_info in enumerate(branch_node.vessel_info):
            segment_vess_ids.append(branch_node.ids[idx])
            segment_gids.append(vess_info['gid'][1]) # out gids
        branch_nodes.append(branch_node)
    
 
        
    
    assert len(junction_gids) == len(junction_outlet_vessels), "Number of junction ids on 3D data will not match the number of outlet vessels in the 0D"
        
    # points of interest
    poi_3d = junction_gids + segment_gids + [0] # include inlet
    poi_0d = lambda df: pd.concat([df.iloc[junction_outlet_vessels]['pressure_in'], df.iloc[segment_vess_ids]['pressure_out'], df.iloc[[0]]['pressure_in']]).to_numpy()
    
    # extract target pressures.
    target_pressures = threed_c.get_pointdata_array("avg_pressure")[poi_3d]
    
    # compute initial case
    tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=True, debug = False)
    init_sim = tmp.run_sim()
    init_sim.convert_to_mmHg()
    pressures_init = poi_0d(init_sim.result_df)

    jcs = []
        
    # submit jobs
    futures = []
    with ProcessPoolExecutor() as executor:
        # iterate through each junction outlet (submit futures)
        for junc_node in junction_nodes:
            for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
                print(f"Changing junction {junc_node.id} downstream vessel {idx}.")
                futures.append(executor.submit(junc_sim, zerod_lpn.get_fast_lpn(), vess, junc_node.id, idx, [junction_outlet_vessels, segment_vess_ids, [0]]))
        # iterate through vessel segments
        for vid in segment_vess_ids:
            print(f"Changing vessel {vid}.")
            futures.append(executor.submit(vess_sim, zerod_lpn.get_fast_lpn(), vid, [junction_outlet_vessels, segment_vess_ids, [0]]))
            
        
        pressures = []
        # parse futures in order
        print("Retrieving results...")
        for idx, f in enumerate(futures):
            
            dr, ps = f.result()
            # compute difference
            pressures.append( ps - pressures_init)
            jcs.append(dr)
            print(f"\tRetrieved results for process {idx}/{len(futures)}")
            
    # convert to numpy
    # add constant & transpose
    pressures.append(list(np.ones(len(pressures[0]))))
    pressures = np.array(pressures).T
    
    # solve for a
    press_inv = np.linalg.inv(pressures) 
    
    # get target pressure differences
    target_pressures_diff = target_pressures - pressures_init
    
    # compute matmul of inverse
    aT = press_inv @ target_pressures_diff
    
    print(target_pressures_diff)
    print(aT)
    #np.save("transform.dat", {'aT': aT, 'target_pressures': target_pressures, 'pressures_init': pressures_init, 'pressures': pressures}, allow_pickle=True)
    
    # iterate through each junction outlet and fill it with appropriate junction values
    counter = 0
    for junc_node in junction_nodes:
        for idx in range(len(junc_node.vessel_info[0]['outlet_vessels'])):
            if aT[counter] > 0: # physical
                zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, R = aT[counter] * jcs[counter], mode = 'add')
            counter += 1
    for vid in segment_vess_ids:
        if aT[counter] > -10:
            zerod_lpn.change_vessel(vessel_id_or_name=vid, R = aT[counter] * jcs[counter], mode = 'add')
        counter += 1
            
    # Split Constant according to Murrays law into proximal resistance
    
    def load_area_file(area_filepath):
        ''' loads a capinfo file
        '''
        with open(area_filepath, 'r') as afile:
            areas = {}
            afile.readline() # ignore first comment line
            for line in afile:
                line = line.rstrip().split()
                areas[line[0]] = float(line[1])
        return areas
    
    def Rpi(Ai, A, Rp):
        return (A / Ai) * Rp
    
    areas= load_area_file(M['workspace']['capinfo'])
    del areas[M['metadata']['inlet']]
    A = sum(list(areas.values()))
    
    # add resistances
    global_const = aT[-1] * m2d(1) / zerod_lpn.inflow.mean_inflow
    for name, bc in zerod_lpn.bc_data.items():
        add_r = Rpi(areas[bc['face_name']], A, global_const)
        print(add_r)
        bc['bc_values']['Rp'] += add_r
    
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
    # reloads the rcrts from previously tuned
    rcr = RCR()
    rcr.read_rcrt_file(M['workspace']['rcrt_file'])
    zerod_lpn.update_rcrs(rcr)
    # write to disk
    zerod_lpn.update()
    
    # load centerlines
    threed_c = Centerlines.load_centerlines(threed_file)
    
    linear_transform(zerod_lpn,threed_c, M)