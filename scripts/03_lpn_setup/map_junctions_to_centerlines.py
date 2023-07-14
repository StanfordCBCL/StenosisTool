# File: map_junctions_to_centerlines.py
# File Created: Monday, 13th February 2023 2:16:14 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 13th July 2023 2:16:04 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Maps the global ids of the centerlines to relevant vessels/junctions on the LPN. Adds Vessels_0D, Junctions_0D, Caps_0D to centerlines and gids to LPN


import argparse

from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.lpn import LPN
from svinterface.manager import Manager
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description= "Maps the global ids of the centerlines to relevant vessels/junctions on the LPN.")
    parser.add_argument("-i", dest = 'config', help = "config.yaml file")
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    # load lpn/centerlines
    lpn = LPN.from_file(M['workspace']['lpn'])
    centerlines = Centerlines.load_centerlines(M['workspace']['centerlines'])
    
    # find GIDS
    centerlines = lpn.find_gids(centerlines)
    
    # rewrite original locations
    centerlines.write_polydata(M['workspace']['centerlines'])
    lpn.update()
    
    
    