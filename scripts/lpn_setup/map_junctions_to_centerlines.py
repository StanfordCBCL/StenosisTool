# File: map_junctions_to_centerlines.py
# File Created: Monday, 13th February 2023 2:16:14 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 14th February 2023 12:39:50 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Finds junctions on centerlines



import argparse
import numpy as np
from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.lpn import LPN
from svinterface.manager import Manager
        
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description= "Adds a Junction array to centerlines.")
    parser.add_argument("-i", dest = 'config', help = "config.yaml file")
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    lpn = LPN.from_file(M['workspace']['lpn'])
    centerlines = Centerlines.load_centerlines(M['workspace']['centerlines'])
    centerlines = lpn.find_gids(centerlines)
    
    # rewrite original location
    centerlines.write_polydata(M['workspace']['centerlines'])
    lpn.update()
    
    
    