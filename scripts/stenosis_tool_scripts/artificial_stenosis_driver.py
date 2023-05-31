# File: artificial_stenosis_driver.py
# File Created: Wednesday, 2nd November 2022 10:30:31 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 3rd April 2023 1:52:08 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Select vessels of a model to generate appropriate amount of stenosis. Generates a stenosis parametrization file containing vessel_id, final r




import argparse
from svinterface.manager.baseManager import Manager
from svinterface.core.zerod.lpn import LPN
from svinterface.utils.io import write_json
from pathlib import Path
from collections import defaultdict
import numpy as np

def random_parametrize(cfg, lpn: LPN):
    
    if not cfg['random']:
        print("Random was not set to true, exitting.")
        exit(1)
    
    # group by generations
    generations = lpn.group_tree_by_generation()
    
    vessels = []
    occlusions = []
    # for every generation specified
    for g, params in cfg['params'].items():
        # get the relevant list of nodes & sort it by increasing length
        nodes = generations[g]
        nodes = sorted(nodes, key = lambda x: x.length(), reverse = True) # sort by greatest length to least
        
        # get total length
        total_len = 0
        for n in nodes:
            total_len+=n.length()
            
        # while the pct is not met, continue occluding
        pct = params[0]
        occluded_len=0
        i = 0
        while occluded_len/total_len < pct:
            vessels.append(nodes[i].ids)
            
            # randomly generate occlusion
            occlusion = np.random.uniform(params[1], params[2])
            for vid in nodes[i].ids:
                lpn.occlude_vessel(vid, occlusion)
            occlusions.append(occlusion)
            
            occluded_len += nodes[i].length()
            i += 1
    return vessels, occlusions
        
         

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Generates artificial stenosis files')
    
    # dev params
    parser.add_argument('-i', dest = 'config', required = True, help = 'config file')
    parser.add_argument('--f', dest = 'force', action = 'store_true', default = False, help = 'Save a new version of artifical stenosis')

    args = parser.parse_args()
    
   
    M = Manager(args.config)
    
    if M['metadata']['diseased'] == True:
        print(" The model is already diseased. Don't add more stenosis to it, regardless of force.")
        exit()
        
    
    print("Generating Artificial Stenosis...", end = '\t', flush = True)
    
    # lpn files
    lpn_dir = Path(M['workspace']['lpn_dir'])
    lpn_file = Path(M['workspace']['lpn'])
    lpn = LPN.from_file(str(lpn_file))
    
    if 'as_simulations' not in M.yaml or M['as_simulations'] is None:
        M.register('as_simulations', {})
    
    
    counter = 0
    # get dir to save it in
    while True:
        rez_dir = lpn_file.parent / (lpn_file.stem + f'.as.{counter}')
        if rez_dir.exists():
            counter += 1
            pass
        else:
            rez_dir.mkdir()
            break
        
    # register dir and files
    M.register(key = counter, value = {}, depth = ['as_simulations'])
    M.register(key = "dir", value = str(rez_dir), depth = ['as_simulations', counter])
    
    # new file
    as_file = rez_dir / (lpn_file.stem + f'.as.in')
    
    ## Get Random Config
    vessel, occ = random_parametrize(M['artificial_stenosis'], lpn)
    
    # convert occlusions to stenosis parametrizations
    sten_p = []
    for o in occ:
        sten_p.append((1/(1-o)) ** (1/2))
    
    
    write_json(str(rez_dir / 'stenosis_parametrization.json'), list(zip(vessel, sten_p)))
    M.register(key="sten_p", value = str(rez_dir / 'stenosis_parametrization.json'), depth = ['as_simulations', counter])

    # write new lpn
    lpn.write_lpn_file(str(as_file))
    M.register(key="lpn", value = str(as_file), depth = ['as_simulations', counter])
    M.register(key="name", value =f"AS-{counter}", depth = ['as_simulations', counter])
        

    M.update()
    print("Done")
    
    
    
    
    
    