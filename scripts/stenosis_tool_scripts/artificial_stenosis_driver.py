# File: artificial_stenosis_driver.py
# File Created: Wednesday, 2nd November 2022 10:30:31 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 27th March 2023 7:12:57 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Select vessels of a model to generate appropriate amount of stenosis. Generates a stenosis parametrization file containing vessel_id, final r




from sgt.utils.parser import ToolParser
from sgt.core.manager import StenosisParametrizationManager
from sgt.core.lpn import LPN
from sgt.utils.io import read_json, write_json

import numpy as np


def use_random(cfg, lpn: LPN):
    ''' randomly generate stenosis
    '''
    gens = cfg['generations']
    # sort by gen number
    gens = sorted(gens, key = lambda x: x[0])
    used_gens = set()
    branch_tree = lpn.get_branch_tree()
    
    vess = []
    occ = []
    for gen, pct_len, occ_low, occ_high in gens:
        # ignore dupe gens
        if gen in used_gens:
            continue
        

        # extract all possible vessels for a gen
        gen_vessels = []
        gen_len = 0
        for node in lpn.tree_bfs_iterator(branch_tree):
            if node.generation == gen:
                gen_vessels.append(node)
                gen_len += node.get_branch_len()

        # randomize which vessels
        indices = np.array(range(len(gen_vessels)))
        np.random.shuffle(indices)
        
        # add vessels to stenosis vessels
        sten_len = 0
        for idx in indices:
            # exit once stenosis length exceeds expected generational pct len
            if sten_len >= gen_len * pct_len:
                break
            n = gen_vessels[idx]
            
            vess.append(n.vess_id)
            occ.append(np.random.uniform(occ_low, occ_high))
            
            sten_len += n.get_branch_len()
            
        
        
        used_gens.add(gen)
    
    return vess, occ

def process_as_config(cfg, lpn: LPN):
    ''' process the vessels specified in config
    '''
    config = read_json(cfg)
    
    # if random, use random
    if config['random'] == True:
        return use_random(config, lpn)
    
    # not random.
    stenosis_vessels = config['vessels']
    vessels = []
    occ = []
    used_vess = set()
    for v in stenosis_vessels:
        # ignore dupes
        if vess_id in used_vess:
            continue
        
        vess_id = v[0]
        occ_low = v[1]
        occ_upper = occ_low
        if len(v) == 3: # if a range is specified
            occ_upper = v[2]

        vessels.append([vess_id])
        occ.append(np.random.uniform(occ_low, occ_upper))
        used_vess.add(vess_id)
    return vessels, occ
        
         

if __name__ == '__main__':
    
    parser = ToolParser(desc = 'Generates artificial stenosis files')
    
    # dev params
    parser.parser.add_argument('-i', help = 'Artificial stenosis config describing which vessels and to what degree to occlude')
    parser.parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'Save a new version of artifical stenosis')

    args = parser.parse_args()
    
   
    SPM = StenosisParametrizationManager(args.config)
    
    if SPM.diseased == True:
        print(" The model is already diseased. Don't add more stenosis to it, regardless of force.")
        exit()
        
    if SPM.SP_lpn.exists():
        if args.force:
            print("Deleting previous artificial stenosis.")
            # unlink previous
            SPM.SP_lpn.unlink(missing_ok = True)
            SPM.stenosis_parametrization.unlink(missing_ok = True )
            SPM.SP_config.unlink(missing_ok = True)
        else:
            print("A previous artificial stenosis model already exists. Doing nothing.")
            exit()
    
    print("Generating Artificial Stenosis...", end = '\t', flush = True)
    lpn = LPN.from_file(str(SPM.lpn))
    
    # proccess config file to obtain vessels
    vessels, occlusion = process_as_config(args.i, lpn)
    
    # modify lpn
    sten_param = {}
    cfg = {'vessels': [], 'random': False}
    for idx, vess_ids in enumerate(vessels):
        sten_param[idx] = [vess_ids, (1/(1-occlusion[idx])) ** (1/2)]
        for vidx in vess_ids:
            lpn.occlude_vessel(vidx, occlusion[idx])
            cfg['vessels'].append([vidx, occlusion[idx]])
        
    # write new lpn
    lpn.write_lpn_file(SPM.SP_lpn)
    
    # write changes_file
    write_json(SPM.stenosis_parametrization, sten_param)
    
    # write config file.
    write_json(SPM.SP_config, cfg, sort_keys = False)
    print("Done")
    
    
    
    
    
    