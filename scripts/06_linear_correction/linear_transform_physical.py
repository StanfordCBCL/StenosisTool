# File: linear_transform_physical.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 17th July 2023 5:12:48 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Perform a linear transform on the junctions, but only saves physical values





import argparse

from svinterface.core.zerod.lpn import LPN, OriginalLPN
from svinterface.core.polydata import Centerlines
from svinterface.core.bc import RCR 
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
from svinterface.manager.baseManager import Manager
from svinterface.utils.misc import m2d
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait



def conc_sim(lpn: OriginalLPN, junc_id: int, which: int, junction_outlet_vessels: list):
    ''' To run simulations using multi-processing '''
    # compute a R as 10%
    r = 1 + lpn.get_junction(junc_id)['junction_values']['R_poiseuille'][which]
    # change to JC
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = r)
    # Solver
    tmp = Solver0Dcpp(lpn, last_cycle_only=True, mean_only=True, debug = False)
    tmp_sim = tmp.run_sim()

    # get results
    tmp_sim.convert_to_mmHg()
    pressures_cur = tmp_sim.result_df.iloc[junction_outlet_vessels]['pressure_in'].to_numpy()
    # undo change
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = 0)

    return pressures_cur       
        
def linear_transform(zerod_lpn: LPN, threed_c: Centerlines, M: Manager):

    # get relevant positions
    tree = zerod_lpn.get_tree()
    # collect
    junction_outlet_vessels = []
    junction_gids = []
    junction_nodes = []
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow='junction'):
        junction_outlet_vessels += junc_node.vessel_info[0]['outlet_vessels']
        junction_gids += junc_node.vessel_info[0]['gid'][1] # out gids
        junction_nodes.append(junc_node)
 
        
    
    assert len(junction_gids) == len(junction_outlet_vessels), "Number of junction ids on 3D data will not match the number of outlet vessels in the 0D"
        
    
    # extract target pressures.
    target_pressures = threed_c.get_pointdata_array("avg_pressure")[junction_gids + [0]]
    
    # compute initial case
    tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=True, debug = False)
    init_sim = tmp.run_sim()
    init_sim.convert_to_mmHg()
    pressures_init = init_sim.result_df.iloc[junction_outlet_vessels + [0]]['pressure_in'].to_numpy()
    
    

        
    # iterate through each junction outlet (submit futures)
    futures = []
    with ProcessPoolExecutor() as executor:
        for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
            for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
                print(f"Changing junction {junc_node.id} vessel {idx}.")
                futures.append(executor.submit(conc_sim, zerod_lpn.get_fast_lpn(), junc_node.id, idx, junction_outlet_vessels + [0]))
        pressures = []
        # parse futures in order
        print("Retrieving results...")
        for idx, f in enumerate(futures):
            
            ps = f.result()
            pressures.append( ps - pressures_init)
            print(f"\tRetrieved results for process {idx}/{len(futures)}")
            
    # convert to numpy
    # add constant & transpose
    pressures.append(list(np.ones(len(pressures[0]))))
    pressures = np.array(pressures).T
    
    # solve for a
    press_inv = np.linalg.inv(pressures) 
    
    # get target pressure differences
    target_pressures_diff = target_pressures - pressures_init
    
    aT = press_inv @ target_pressures_diff
    
    
    print(target_pressures_diff)
    print(aT)
    # np.save("transform.dat", {'aT': aT, 'target_pressures': target_pressures, 'pressures_init': pressures_init, 'pressures': pressures}, allow_pickle=True)
    
    # iterate through each junction outlet and fill it with appropriate junction values
    counter = 0
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
        for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
            v = zerod_lpn.get_vessel(vess)
            if aT[counter] + junc_node.vessel_info[0]['junction_values']['R_poiseuille'][idx] > -v['zero_d_element_values']['R_poiseuille']: # physical
                zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, R = aT[counter], mode = 'add')
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
    
    # add resistances in equal ratio
    global_const = aT[-1] * m2d(1) / zerod_lpn.inflow.mean_inflow
    for name, bc in zerod_lpn.bc_data.items():
        add_r = Rpi(areas[bc['face_name']], A, global_const)
        print(add_r)
        rat = bc['bc_values']['Rp'] / (bc['bc_values']['Rd'] + bc['bc_values']['Rp'])
        bc['bc_values']['Rp'] += add_r * rat
        bc['bc_values']['Rd'] += add_r * (1-rat)
    
    # save the lpn.
    zerod_lpn.update()
    
    

    

    
    
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Perform a linear optimization on the branches")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    
    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    zerod_file = M['workspace']['base_lpn']
    threed_file = M['workspace']['3D']
    
    # reset lpn back to base
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.write_lpn_file(M['workspace']['lpn'])
    # re-init LPN
    zerod_lpn = LPN.from_file(M['workspace']['lpn'])

    
    # load centerlines
    threed_c = Centerlines.load_centerlines(threed_file)
    
    linear_transform(zerod_lpn,threed_c, M)