# File: stenosis_detection_driver.py
# File Created: Wednesday, 29th June 2022 11:26:57 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 4th April 2023 11:34:37 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Performs test of pressure drop across tree


from svinterface.core.zerod.lpn import LPN
from svinterface.core.zerod.solver import Solver0Dcpp, SolverResults
from svinterface.utils.io import write_json
from svinterface.manager.baseManager import Manager
import numpy as np
import argparse



def detect_stenosis(lpn: LPN, threshold ):
    ''' determines stenosis parametrization
    '''
    
    sten_param = {}
    solver = Solver0Dcpp(lpn, last_cycle_only=True, debug = False)
    solver_results = solver.run_sim()

    counter = 0
    # compute pressure drops.
    mpa = lpn.get_tree()
    for node in lpn.tree_bfs_iterator(mpa, "branch"):
        
        # check the vessel
        pin = solver_results.get_avg_val(node.vessel_info[0]['vessel_name'], 'pressure_in')
        pout = solver_results.get_avg_val(node.vessel_info[-1]['vessel_name'], 'pressure_out')
        print(" Node id:", node.ids[0], "Pin/Pout:", pin, pout,"Ratio:", pout/pin)
        if pout/pin < (1 - threshold):
            repair_area_increase = (1/((((threshold)) / (1 - pout/pin)) ** 1/4))
            print("Repair:", repair_area_increase)
            sten_param[counter] = [node.ids, repair_area_increase]
            # modify lpn
            for i in node.ids:
                lpn.repair_vessel(i, repair_area_increase)
            counter += 1
    return sten_param

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Detects stenosis in a diseased model')
    
    # dev params
    parser.add_argument('-i', dest = 'config', help = 'Yaml File')
    parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force a new iteration over the old')
    parser.add_argument('-t', type = float, default = .2,  help = 'determine the pressure drop threshold that counts as stenosis: default = .1 (10%)')
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    
    print("Detecting Stenosis regions...", end = '\t', flush = True)
    lpn = LPN.from_file(M['workspace']['lpn'])

    sten_param = detect_stenosis(lpn, args.t)
    
    