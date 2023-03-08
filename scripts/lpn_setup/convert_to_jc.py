# File: convert_to_jc.py
# File Created: Tuesday, 1st November 2022 1:28:49 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 27th February 2023 11:33:02 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Postprocess a 0D model to add Junction coefficients.

from svinterface.core.zerod.lpn import LPN
from svinterface.manager.baseManager import Manager
import argparse

def convert_to_jc(lpn: LPN):
    
    for junc in lpn.junctions:
        if junc['junction_type'] != 'internal_junction':
            
            S0 = junc['areas'][0]
            s_outlets = junc['areas'][1:]

            
            density = lpn.simulation_params['density']

            # update
            for idx, S1 in enumerate(s_outlets):
                jc = 1.9380249380494794 * density * ( (S0/S1 - 1) **2) / (2 * S0**2)
                lpn.change_junction_outlet(junc['junction_name'], which = idx, S = jc)
                
           

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=  "Convert a normal LPN to a JC model version.")
    parser.add_argument("-i", dest = 'config', help = "config.yaml file")
    
    args = parser.parse_args()
    
    M = Manager(args.config)

    print('Converting to JC LPN...', end = '\t', flush = True)
 
    # conversion
    lpn_jc = LPN.from_file(str(M['workspace']['lpn']))
    lpn_jc.to_cpp(normal = False)
    
    convert_to_jc(lpn_jc)
    
    # write new file
    lpn_jc.update()
    
    ## Cleanup

    print('Done')