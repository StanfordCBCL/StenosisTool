# File: map_0D_to_centerlines.py
# File Created: Monday, 23rd January 2023 7:16:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 3rd April 2023 2:14:17 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Takes 0D results csv file and corresponding centerlines used to construct 0D LPN, maps a 0D model's results to centerlines.

from svinterface.core.zerod.lpn import LPN
from svinterface.core.polydata import Polydata
from svinterface.core.zerod.solver import SolverResults
from svinterface.manager.baseManager import Manager
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
    parser.add_argument("-i", dest = 'config', required = True, help = 'config file')
    
    parser.add_argument('-mode', default = None, help = 'Mode: None, AS, R')
    parser.add_argument("-sim", type = int, help = "Simulation number to use")

    # flags
    parser.add_argument("--mmHg", action = 'store_true', default = False, help ='converts all values to mmHg')
    parser.add_argument("--s", dest = 'summary', action = 'store_true', default = False, help ='only reports summary values only')

    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    
    # load original centerlines
    c = Polydata.load_polydata(M['workspace']['centerlines'])
    
    sims = 'simulations'
    if args.mode == 'AS':
        sims = 'as_simulations'
    elif args.mode == 'R':
        sims = 'r_simulations'
    else:
        raise ValueError("-mode must be AS or R or not set")
        
    # load LPN
    lpn = LPN.from_file(M[sims][args.sim]['lpn'])
    
    # load results
    results = SolverResults.from_csv(M[sims][args.sim]['csv'])
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
    out_poly = str(Path(M[sims][args.sim]['dir']) / output_file)
    c.write_polydata(out_poly)
    
    M.register(key = 'centerlines', value = out_poly, depth = [sims, args.sim])
    M.update()