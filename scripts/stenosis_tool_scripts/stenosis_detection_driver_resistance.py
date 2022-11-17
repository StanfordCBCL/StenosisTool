# File: stenosis_detection_driver.py
# File Created: Wednesday, 29th June 2022 11:26:57 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 3rd November 2022 12:39:56 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Performs test of pressure drop across tree


from sgt.core.manager import StenosisParametrizationManager
from sgt.utils.parser import ToolParser
from sgt.core.lpn import LPN
from sgt.utils.io import write_json
import numpy as np

def diameter_formula(order, age):
    ''' from https://journals.physiology.org/doi/full/10.1152/ajpheart.00123.2020 
    Equation for diameter '''
    a = 1.203e-4
    b = 0.3927
    c = 0.2381
    d = 0.001 + a * order * np.exp(b * order) * (age ** c)
    return d

def detect_stenosis(lpn: LPN, occlusion, max_gens, age ):
    ''' determines stenosis parametrization
    '''
    tree = lpn.get_branch_tree()
    
    sten_param = {}
    counter = 0
    for node in lpn.tree_bfs_iterator(tree):
        if node.generation < max_gens:
            control_radius = diameter_formula(16 - node.generation, age) / 2
           
            for vess_idx in node.vess_id:
                true_radius = lpn.get_vessel_radius(vess_idx)
                print(node.generation ,control_radius, true_radius, 1 - (true_radius / control_radius) ** 2,occlusion)
                #print(1 - (true_radius / control_radius) ** 2)
                if (1 - (true_radius / control_radius) ** 2) > occlusion:
                    sten_param[counter] = [[node.vess_id], control_radius]
                    counter += 1
    print(sten_param)
        
    return sten_param

if __name__ == '__main__':
    
    parser = ToolParser(desc='Detects stenosis in a diseased model')
    
    # dev params
    parser.parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force a new iteration over the old')
    parser.parser.add_argument('-occ', type = float, default = .7,  help = 'determine the occlusion that counts as stenosis')
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
    lpn.write_lpn_file(SPM.SP_lpn)
    
    sten_param = detect_stenosis(lpn, args.occ, args.g, SPM.info['metadata']['age'])
    
    write_json(SPM.stenosis_parametrization, sten_param)
    print("Done")
    
    