# File: convert_to_jc.py
# File Created: Tuesday, 1st November 2022 1:28:49 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 6th February 2023 12:12:09 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Postprocess a 0D model to add Junction coefficients.

from sgt.core.lpn import LPN
from sgt.core.manager import LPNConstructionManager
from sgt.utils.parser import ToolParser


def convert_to_jc(lpn: LPN):
    
    for junc in lpn.junctions:
        if junc['junction_type'] != 'internal_junction':
            
            inlet = junc['inlet_vessels']
            outlets = junc['outlet_vessels']
            S0 = junc['areas'][0]
            s_outlets = junc['areas'][1:]
            
            assert len(outlets) == len(s_outlets), 'Number of areas does not match number of outlets'
            
            outlet_areas = list(zip(outlets, s_outlets))
            
            density = lpn.simulation_params['density']
            # rename
            junc['junction_type'] = "BloodVesselJunction"
            # update
            junc['junction_values'] = {"stenosis_coefficient": []}
            for out, S1 in outlet_areas:
                jc = 1.52 * density * ( (S0/S1 - 1) **2) / (2 * S0**2)
                junc['junction_values']['stenosis_coefficient'].append(jc)
           

if __name__ == '__main__':
    parser = ToolParser(desc = "Convert a normal LPN to a JC model version.")
    args = parser.parse_args()
    
    M = LPNConstructionManager(args.config)
    
    if M.jc:
        print("Already a JC model. No action was taken.")
        exit(0)

    print('Converting to JC LPN...', end = '\t', flush = True)
 
    # conversion
    lpn_jc = LPN.from_file(str(M.lpn))
    convert_to_jc(lpn_jc)
    
    # add jc to config
    M.config_add(['options','jc'], True)
    M.write_config()
    
    # delete old file
    M.lpn.unlink()
    
    # reinit
    M = LPNConstructionManager(args.config)
    
    # write new file
    lpn_jc.write_lpn_file(M.lpn)
    
    ## Cleanup

    print('Done')