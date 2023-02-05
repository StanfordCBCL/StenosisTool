# File: map_0D_to_centerlines.py
# File Created: Monday, 23rd January 2023 7:16:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 26th January 2023 5:54:20 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Takes 0D results csv file and corresponding centerlines used to construct 0D LPN, maps a 0D model's results to centerlines.

from sgt.core.lpn import LPN
from sgt.core.polydata import Polydata
from sgt.core.solver import SolverResults
import argparse
from pathlib import Path
import numpy as np

def remove_nonsummary(centerlines: Polydata):
    
    for arr_name in centerlines.get_pointdata_arraynames():
        if arr_name.startswith("pressure") or arr_name.startswith("flow"):
            centerlines.remove_pointdata_array(arr_name)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Takes 0D results csv file and corresponding centerlines used to construct 0D LPN, maps a 0D model's results to centerlines.")
    
    # files
    parser.add_argument("-c", dest = 'centerlines', required = True, help = 'Centerlines file')
    parser.add_argument("-r", dest = "results",required = True, help = 'results csv file')
    parser.add_argument("-lpn", required = True, help = "Lumped Param Network")
    parser.add_argument("-o", required = True, help = 'output file of centerlines')
    # flags
    parser.add_argument("--mmHg", action = 'store_true', default = False, help ='converts all values to mmHg')
    parser.add_argument("--s", dest = 'summary', action = 'store_true', default = False, help ='only reports summary values only')

    args = parser.parse_args()
    
    # load original centerlines
    c = Polydata()
    c.load_polydata(args.centerlines)
    
    # load LPN
    lpn = LPN.from_file(args.lpn)
    
    # load results
    results = SolverResults.from_csv(args.results)
    if args.mmHg:
        results.convert_to_mmHg()

    # update centerlines
    c = results.project_to_centerline(lpn, c)
    
    output_file = 'centerline_projection.vtp'
    
    # keep summary data only
    if args.summary:
        remove_nonsummary(c)
        output_file = 'centerline_projection.summary.vtp'
        
    # write polydata
    c.write_polydata(str(Path(args.o) / output_file))