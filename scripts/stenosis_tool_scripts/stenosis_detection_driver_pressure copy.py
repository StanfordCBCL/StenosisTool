# File: stenosis_detection_driver.py
# File Created: Wednesday, 29th June 2022 11:26:57 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 3rd April 2023 11:04:48 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Performs test of pressure drop across tree


from 
from sgt.utils.parser import ToolParser
from sgt.core.lpn import LPN
from sgt.core.solver import Solver0Dcpp, SolverResults
from sgt.utils.io import write_json
import numpy as np



def detect_stenosis(lpn: LPN, threshold, max_gens ):
    ''' determines stenosis parametrization
    '''
    
    sten_param = {}
    solver = Solver0Dcpp(lpn, last_cycle_only=True, debug = False)
    solver_results = solver.run_sim()

    counter = 0
    # compute pressure drops.
    mpa = lpn.get_vessel_tree()
    for node in lpn.tree_bfs_iterator(mpa):
        # skip if generation is too high
        if node.generation > max_gens:
            continue
        
        # check the vessel
        pin = solver_results.get_avg_val(node.vessel_info[0]['vessel_name'], 'pressure_in')
        pout = solver_results.get_avg_val(node.vessel_info[0]['vessel_name'], 'pressure_out')
        print("Generation: ", node.generation, "Node id:", node.vess_id[0], "Pin/Pout:", pin, pout,"Ratio:", pout/pin)
        if pout/pin < (1 - threshold):
            repair_area_increase = (1/((((threshold)) / (1 - pout/pin)) ** 1/4))
            print("Repair:", repair_area_increase)
            sten_param[counter] = [node.vess_id, repair_area_increase]
            # modify lpn
            lpn.repair_vessel(node.vess_id[0], repair_area_increase)
            counter += 1
    return sten_param

if __name__ == '__main__':
    
    parser = ToolParser(desc='Detects stenosis in a diseased model')
    
    # dev params
    parser.parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force a new iteration over the old')
    parser.parser.add_argument('-t', type = float, default = .2,  help = 'determine the pressure drop threshold that counts as stenosis: default = .1 (10%)')
    parser.parser.add_argument('-g', type = int,  default = 3, help = 'The max number of generations to use.')
    args = parser.parse_args()
    
    SPM = StenosisParametrizationManager(args.config)
    
    if SPM.diseased == False:
        print("This models is not diseased and does not require any repair.")
        exit(1)
    
    if SPM.SP_lpn.exists():
        if args.force:
            print("Deleting previous stenosis detection.")
            # unlink previous
            SPM.SP_lpn.unlink(missing_ok = True)
            SPM.stenosis_parametrization.unlink(missing_ok = True )
        else:
            print("Stenosis detection was already performed. Nothing occured.")
            exit(1)
    
    print("Detecting Stenosis regions...", end = '\t', flush = True)
    lpn = LPN.from_file(SPM.lpn)

    sten_param = detect_stenosis(lpn, args.t, args.g)
    lpn.write_lpn_file(SPM.SP_lpn)
    
    write_json(SPM.stenosis_parametrization, sten_param)
    print("Done")
    
    