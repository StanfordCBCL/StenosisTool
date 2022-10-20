# File: add_stenosis_locations.py
# File Created: Wednesday, 27th July 2022 11:28:30 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 18th October 2022 9:14:13 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Adds stenosis locations to the centerline files

import argparse
from src.polydata import Centerlines
from src.file_io import read_json
from src.misc import create_tool_parser, get_solver_path
from src.lpn import Solver0D
import numpy as np
import re
from pathlib import Path

def main(args):
    
    for solver_dir in args.solver_dirs:
        
        solver_dir = Path(solver_dir)
        
        
        # centerlines
        centerlines = Centerlines()
        centerlines.load_centerlines(solver_dir / 'model_centerlines.vtp')
        # stenosis file
        stenosis_info = read_json(solver_dir / 'stenosis_vessels.dat')
        # solver file
        solver_file = get_solver_path(solver_dir)
        solver0d = Solver0D()
        solver0d.read_solver_file(solver_file)        
        
        # Number stenosis
        stenosis_vessels = {idx:vess for idx, vess in enumerate(stenosis_info['all_changed_vessels'])}
        
        # retrieve branch ids and create stenosis array
        branch_ids = centerlines.get_pointdata(Centerlines.PointDataFields.BRANCHID)
        stenosis_nodes = np.ones_like(branch_ids) * -1

        for (idx, vessels) in stenosis_vessels.items():
            # find branch id of vessels
            v = solver0d.get_vessel(vessels[0])
            branch_id = re.search("branch([0-9]+)_.*", v['vessel_name']).group(1)
            stenosis_nodes = np.where(branch_ids == int(branch_id), idx, stenosis_nodes )
        
        centerlines.add_pointdata(stenosis_nodes, 'Stenosis Vessels')
        centerlines.write_centerlines(solver_dir / 'model_centerlines.vtp')

if __name__ == '__main__':
    
    parser = create_tool_parser('Add stenosis locations')
    args = parser.parse_args()
    
    main(args)