

import argparse
import numpy as np
from sgt.core.polydata import Centerlines
from sgt.core.lpn import LPN
        
        

if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description= "Adds a Junction array to centerlines.")
    
    parser.add_argument("-c", dest = 'centerlines', help = 'centerline file')
    parser.add_argument("-lpn", help = "LPN")
    
    args = parser.parse_args()
    
    lpn = LPN.from_file(args.lpn)
    centerlines = Centerlines.from_file(args.centerlines)
    centerlines = lpn.find_junctions_on_centerlines(centerlines)
    
    # rewrite original location
    centerlines.write_polydata(args.centerlines)
    lpn.write_lpn_file(lpn.lpn_file)
    
    
    