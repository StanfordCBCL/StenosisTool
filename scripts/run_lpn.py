# File: run_lpn.py
# File Created: Thursday, 3rd November 2022 12:44:51 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 3rd November 2022 1:43:42 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Solves an LPN

from sgt.core.solver import Solver0Dcpp
from sgt.core.lpn import LPN
from sgt.utils.parser import Parser

from pathlib import Path

if __name__ == '__main__':
    
    parser = Parser(desc = 'Solves an LPN')
    
    parser.parser.add_argument('-i', dest = 'lpn', help = 'path to input lpn')
    parser.parser.add_argument('-c', dest = 'csv', action = 'store_true', default = False, help = 'save csv file: Default = False')
    parser.parser.add_argument('-b', dest = 'branch', action = 'store_true', default=False,  help = 'to convert c output to python branch files: Default = False')
    parser.parser.add_argument('-l', dest = 'last_cycle', action = 'store_true', default = False, help = 'only save the last cycle worth of results: Default = False')
    parser.parser.add_argument('-m', dest = 'mean_only', action = 'store_true', default = False, help = 'only save the mean of the results: Default = False')
    parser.parser.add_argument('-v', dest = 'validate', action = 'store_true', default = False, help = 'validate the run with inlet pressure waveform: Default = False')
    
    
    args = parser.parse_args()
    
    lpn_file =  Path(args.lpn)
    if not lpn_file.exists():
        raise FileNotFoundError(f'LPN File {str(lpn_file)} not found.')
    
    lpn = LPN.from_file(str(lpn_file))
    
    # get dir to save it in
    counter = 0
    while True:
        rez_dir = lpn_file.parent / (lpn_file.stem + f'.sim.{counter}')
        if rez_dir.exists():
            counter += 1
        else:
            rez_dir.mkdir()
            break
    
    solver = Solver0Dcpp(lpn, last_cycle_only=args.last_cycle, mean_only=args.mean_only, debug = True)
    
    results = solver.run_sim_pipeline(validate = args.validate, save_csv = args.csv, save_branch = args.branch, out_dir = rez_dir)
    
    
    