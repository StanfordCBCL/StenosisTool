# File: linear_transform.py
# File Created: Tuesday, 14th February 2023 11:25:35 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 15th March 2023 3:09:38 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description:  Perform a linear transform on the junctions 
#! USES SYS, DIA, AND AVG






import argparse

from svinterface.core.zerod.lpn import LPN, OriginalLPN
from svinterface.core.polydata import Centerlines
from svinterface.core.bc import RCR 
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
from svinterface.manager.baseManager import Manager
from svinterface.utils.misc import m2d
import numpy as np
from concurrent.futures import ProcessPoolExecutor, wait



def conc_sim(lpn: OriginalLPN, vess: int, junc_id: int, which: int, junction_outlet_vessels: list):
    ''' To run simulations using multi-processing '''
    # get downstream vessel.
    v = lpn.get_vessel(vess)
    # compute a R as 10%
    r = .1 * v['zero_d_element_values']['R_poiseuille']
    # change to JC
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = r)
    # Solver
    tmp = Solver0Dcpp(lpn, last_cycle_only=True, mean_only=False, debug = False)
    tmp_sim = tmp.run_sim()

    # get results
    tmp_sim.convert_to_mmHg()
    
    pressures_cur = [tmp_sim.get_avg_val(lpn.get_vessel(0)['vessel_name'], "pressure_in")] 
    
    p_cur_mean = []
    p_cur_sys = []
    p_cur_dia = []
    for vid in junction_outlet_vessels:
        vess_name = lpn.get_vessel(vid)['vessel_name']
        p_cur_mean.append(tmp_sim.get_avg_val(vess_name, "pressure_in"))
        p_cur_sys.append(tmp_sim.get_max_val(vess_name, "pressure_in"))
        p_cur_dia.append(tmp_sim.get_min_val(vess_name, "pressure_in"))
    pressures_cur = np.array(pressures_cur + p_cur_mean + p_cur_sys + p_cur_dia)
    # undo change
    lpn.change_junction_outlet(junction_id_or_name = junc_id, which = which, R = 0)

    return r, pressures_cur       
        
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
    target_pressures = threed_c.get_pointdata_array("avg_pressure")[[0] + junction_gids]
    for arr_name in threed_c.get_pointdata_arraynames():
        if arr_name.startswith("sys_pressure"):
            sys_name = arr_name
        elif arr_name.startswith("dia_pressure"):
            dia_name = arr_name
    target_pressures = np.append(target_pressures, threed_c.get_pointdata_array(sys_name)[junction_gids])
    target_pressures = np.append(target_pressures, threed_c.get_pointdata_array(dia_name)[junction_gids])
    
    # compute initial case
    tmp = Solver0Dcpp(zerod_lpn, last_cycle_only=True, mean_only=False, debug = False)
    init_sim = tmp.run_sim()
    init_sim.convert_to_mmHg()
    
    # MPA
    pressures_init = [init_sim.get_avg_val(zerod_lpn.get_vessel(0)['vessel_name'], "pressure_in")] 
    
    p_init_mean = []
    p_init_sys = []
    p_init_dia = []
    for vid in junction_outlet_vessels:
        vess_name = zerod_lpn.get_vessel(vid)['vessel_name']
        p_init_mean.append(init_sim.get_avg_val(vess_name, "pressure_in"))
        p_init_sys.append(init_sim.get_max_val(vess_name, "pressure_in"))
        p_init_dia.append(init_sim.get_min_val(vess_name, "pressure_in"))
    pressures_init = np.array(pressures_init + p_init_mean + p_init_sys + p_init_dia)
    
    jcs = []
    # iterate through each junction outlet (submit futures)
    futures = []
    with ProcessPoolExecutor() as executor:
        for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
            for idx, vess in enumerate(junc_node.vessel_info[0]['outlet_vessels']):
                print(f"Changing junction {junc_node.id} vessel {idx}.")
                futures.append(executor.submit(conc_sim, zerod_lpn.get_fast_lpn(), vess, junc_node.id, idx, junction_outlet_vessels))
        pressures = [list(np.ones(len(pressures_init)))]
        # parse futures in order
        print("Retrieving results...")
        for idx, f in enumerate(futures):
            
            r, ps = f.result()
            pressures.append( ps - pressures_init)
            jcs.append(r)
            print(f"\tRetrieved results for process {idx}/{len(futures)}")
            
    # convert to numpy
    # add constant & transpose
    pressures = np.array(pressures).T
    
    print(pressures.shape)
    print(pressures_init.shape)
    print(target_pressures.shape)
    
    # solve for a
    press_inv = np.linalg.pinv(pressures) 
    
    # get target pressure differences
    target_pressures_diff = target_pressures - pressures_init
    
    aT = press_inv @ target_pressures_diff
    
    
    #! SAVING STUFF FOR DEBUGGING
    print(target_pressures_diff)
    print(aT)
    np.save("transform.dat", {'aT': aT, 'target_pressures': target_pressures, 'pressures_init': pressures_init, 'pressures': pressures}, allow_pickle=True)
    
    # iterate through each junction outlet and fill it with appropriate junction values
    counter = 1
    for junc_node in zerod_lpn.tree_bfs_iterator(tree, allow = 'junction'):
        for idx in range(len(junc_node.vessel_info[0]['outlet_vessels'])):
            if aT[counter] > -10: # physical
                zerod_lpn.change_junction_outlet(junction_id_or_name=junc_node.id, which=idx, R = aT[counter] * jcs[counter])
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
    global_const = aT[0] * m2d(1) / zerod_lpn.inflow.mean_inflow
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