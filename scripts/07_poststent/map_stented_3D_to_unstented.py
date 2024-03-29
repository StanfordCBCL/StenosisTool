# File: map_stented3D_to_unstented.py
# File Created: Tuesday, 4th July 2023 12:21:05 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:26:21 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Maps stented 3D extracted centerlines to the diseased centerlines

from svinterface.core.polydata import Centerlines
from svinterface.manager import Manager

import argparse
import re
import numpy as np
from pathlib import Path

def map_stented_to_prestent(stented: Centerlines, prestent: Centerlines):
    """ Performs map of stented to prestented"""
    j = stented.get_pointdata_array('Junctions_0D')
    v = stented.get_pointdata_array('Vessels_0D')
    c = stented.get_pointdata_array('Caps_0D')
    
    def add_valid(p_arr, s_arr):
        """Find valid regions"""
        jidx = np.where(j > -1)
        p_arr[j[jidx]] = s_arr[jidx]
        vidx = np.where(v > -1)
        p_arr[v[vidx]] = s_arr[vidx]
        cidx = np.where(c > -1)
        p_arr[c[cidx]] = s_arr[cidx]
        return p_arr
    
    pgid = prestent.get_pointdata_array("GlobalNodeId")
    
    
    for f in ['pressure', 'velocity']:      
        for arr_name in stented.get_pointdata_arraynames():   
            x = re.search(f"{f}_([0-9]+)", arr_name)
            #! Assumes in time order
            if x:
               p_arr = add_valid(np.zeros_like(pgid), stented.get_pointdata_array(arr_name))
               prestent.add_pointdata(p_arr, array_name=arr_name)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="map stented to prestent")
    parser.add_argument("-c", dest='stent', help = 'stented 3D extracted centerlines (must have been matched via centerline_match.py)' )
    parser.add_argument("-i", dest='config', help= 'Config file for project')
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    # load centerlines
    p = Centerlines.load_centerlines(M['workspace']['centerlines'])
    s = Centerlines.load_centerlines(args.stent)
    
    # map
    map_stented_to_prestent(s, p)
    
    # write mapped centerlines back to same dir
    opath = Path(args.stent)
    p.write_polydata(opath.parent / (str(opath.stem) + ".mapped.vtp" ))
    