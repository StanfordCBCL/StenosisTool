# File: linear_transform_global_split.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:29:10 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description:  Perform a linear correction all resistances in diseased model


from svinterface.core.zerod.lpn import LPN, OriginalLPN
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.solver import Solver0Dcpp
from svinterface.manager.baseManager import Manager
from svinterface.utils.misc import m2d

import argparse
import numpy as np
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import json



def junc_sim(lpn: OriginalLPN, junc_id: int, which: int, poi_0d):
    ''' To run simulations using multi-processing '''
    
    # get downstream vessel.
    r = 1
    # change to JC
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = r, mode='add')
    # Solver
    tmp = Solver0Dcpp(lpn, last_cycle_only=True, mean_only=True, debug = False)
    tmp_sim = tmp.run_sim()

    # get results
    tmp_sim.convert_to_mmHg()
    pressures_cur = pd.concat([tmp_sim.result_df.iloc[poi_0d[0]]['pressure_in'], tmp_sim.result_df.iloc[poi_0d[1]]['pressure_out'], tmp_sim.result_df.iloc[poi_0d[2]]['pressure_in']]).to_numpy()
    # undo change
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = -r, mode = 'add')

    return pressures_cur   

def vess_sim(lpn: OriginalLPN, vess: int, poi_0d):
    ''' To run simulations using multi-processing '''
    dr = 1
    # change to JC
    lpn.change_vessel(vessel_id = vess, R = dr, mode = 'add')
    # Solver
    tmp = Solver0Dcpp(lpn, last_cycle_only=True, mean_only=True, debug = False)
    tmp_sim = tmp.run_sim()

    # get results
    tmp_sim.convert_to_mmHg()
    pressures_cur = pd.concat([tmp_sim.result_df.iloc[poi_0d[0]]['pressure_in'], tmp_sim.result_df.iloc[poi_0d[1]]['pressure_out'],  tmp_sim.result_df.iloc[poi_0d[2]]['pressure_in']]).to_numpy()
    # undo change
    lpn.change_vessel(vessel_id = vess, R = -dr, mode = 'add')

    return pressures_cur       
        
def linear_transform(zerod_lpn: LPN, threed_c: Centerlines, M: Manager,iter: int):
    
    # get relevant positions
    tree = zerod_lpn.get_tree()
    # determine sides
    zerod_lpn.det_lpa_rpa(tree)
    
    total_glob = 0
    for i in range(iter):
        for side in  'MPA', 'RPA', 'LPA':
            print(f"Evaluating {side}.")
            total_glob += linear_transform_side(zerod_lpn, threed_c, M, side, i)
    print(total_glob)
        
def linear_transform_side(zerod_lpn: LPN, threed_c: Centerlines, M: Manager, side, iter):

    # get relevant positions
    tree = zerod_lpn.get_tree()
    # collect junction values
    junction_outlet_vessels = []
    junction_gids = []
    junction_nodes = []
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow='junction'):
        if junc_node.vessel_info[0]['side'] == side:
            junction_outlet_vessels += junc_node.vessel_info[0]['outlet_vessels']
            junction_gids += junc_node.vessel_info[0]['gid'][1] # out gids
            junction_nodes.append(junc_node)
        
    # collect branch values
    segment_vess_ids = []
    segment_gids = []
    branch_nodes = []
    for branch_node in zerod_lpn.tree_bfs_iterator(tree, allow='branch'):
        if branch_node.vessel_info[0]['side'] == side:
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

        
    # submit jobs
    futures = []
    with ProcessPoolExecutor() as executor:
        # iterate through each junction outlet (submit futures)
        for junc_node in junction_nodes:
            for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
                print(f"Changing junction {junc_node.id} downstream vessel {idx}.")
                futures.append(executor.submit(junc_sim, zerod_lpn.get_fast_lpn(), junc_node.id, idx, [junction_outlet_vessels, segment_vess_ids, [0]]))
        # iterate through vessel segments
        for vid in segment_vess_ids:
            print(f"Changing vessel {vid}.")
            futures.append(executor.submit(vess_sim, zerod_lpn.get_fast_lpn(), vid, [junction_outlet_vessels, segment_vess_ids, [0]]))
            
        
        pressures = []
        # parse futures in order
        print("Retrieving results...")
        for idx, f in enumerate(futures):
            
            ps = f.result()
            #print(ps, pressures_init)
            # compute difference
            pressures.append( ps - pressures_init)
            print(f"\tRetrieved results for process {idx}/{len(futures)}")
            
    # convert to numpy
    # add constant & transpose
    pressures.append(np.ones_like(pressures[-1]))
    pressures = np.array(pressures).T
    
    # solve for a
    press_inv = np.linalg.inv(pressures) 

    # get target pressure differences
    target_pressures_diff = target_pressures - pressures_init
    
    # compute matmul of inverse
    aT = press_inv @ target_pressures_diff
    
    # iterate through each junction outlet and fill it with appropriate junction values
    counter = 0
    for junc_node in junction_nodes:
        for idx in range(len(junc_node.vessel_info[0]['outlet_vessels'])):
            if aT[counter] + junc_node.vessel_info[0]['junction_values']['R_poiseuille'][idx] > 0: # physical
                zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, R = aT[counter], mode = 'add')
            else:
                # otherwise minimum which is 0
                zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, R = 0)
            counter += 1
    for vid in segment_vess_ids:
        v = zerod_lpn.get_vessel(vid)
        if aT[counter] > -v['zero_d_element_values']['R_poiseuille']: # physical
            zerod_lpn.change_vessel(vessel_id_or_name=vid, R = aT[counter], mode = 'add')
        else:
            zerod_lpn.change_vessel(vessel_id_or_name=vid, R = 0)
        counter += 1
    
    
    # Split Constant according to Murrays law into proximal resistance and distal resistance
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
    
    # add resistances if the iteration isn't first, otherwise don't add it. (Otherwise, the global constant in the MPA initially may cause negatives in BC's making the resulting network unsolvable)
    if iter != 0:
        global_const = aT[-1] * m2d(1) / zerod_lpn.inflow.mean_inflow
        for name, bc in zerod_lpn.bc_data.items():
            add_r = Rpi(areas[bc['face_name']], A, global_const)
            # print(add_r)
            rat = bc['bc_values']['Rp'] / (bc['bc_values']['Rd'] + bc['bc_values']['Rp'])
            bc['bc_values']['Rp'] += add_r * rat
            bc['bc_values']['Rd'] += add_r * (1-rat)
    
    # save the lpn.
    zerod_lpn.update()
    
    return aT[-1]
    
    

    

    
    
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Perform a second linear optimization on the branches")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    parser.add_argument("-c", dest = 'cent', help = 'mapped formatted centerlines of stented models')
    parser.add_argument("-n", dest = 'name', help = 'Name for stented model')
    parser.add_argument("--iter", dest = 'iter', type=int, default = 5, help = 'Number of iterations to run correction')
    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    # get 0D LPN and 3D cent
    zerod_file = M['parameterization']['base_lpn']
    threed_file = args.cent
    
    # correction dir
    correction_dir = Path(M['workspace']['param_dir']) / args.name
    correction_dir.mkdir(exist_ok = True)
    M.register(args.name, {}, ['parameterization', 'corrections'])
    M.register('dir', str(correction_dir), ['parameterization','corrections', args.name])
    
    # copy lpn to dir
    new_lpn_path = correction_dir / Path(zerod_file).name
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.write_lpn_file(str(new_lpn_path))
    zerod_lpn = LPN.from_file(str(new_lpn_path))
    M.register('lpn', str(new_lpn_path), ['parameterization','corrections', args.name])

    # reinit LPN
    zerod_lpn = LPN.from_file(str(new_lpn_path))
    
    # add that relevant regions is all vessels and junctions
    vess_ids = list(range(zerod_lpn.num_vessels()))
    junc_ids = ['J' + str(i) for i in range(zerod_lpn.num_junctions())]
    relevant_regions = correction_dir / 'relevant_regions.json'
    with relevant_regions.open("w") as rrfile:
        json.dump({'Vessels': sorted(list(vess_ids)), 'Junctions': sorted(list(junc_ids))}, rrfile, indent=4, sort_keys=True)
    M.register('relevant_regions', str(relevant_regions), ['parameterization','corrections', args.name])
    
    # load centerlines
    threed_c = Centerlines.load_centerlines(threed_file)
    
    # linear transform
    linear_transform(zerod_lpn, threed_c, M, iter = args.iter)
    
    # run the pipeline
    solver = Solver0Dcpp(zerod_lpn, debug = True)
    solver.run_sim_pipeline(validate=True, save_csv=True, save_branch=False, out_dir = str(correction_dir))
    M.register('csv', str(correction_dir / 'branch_results.csv'), ['parameterization','corrections', args.name])
    M.update()
    