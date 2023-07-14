# File: run_lpn.py
# File Created: Thursday, 3rd November 2022 12:44:51 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 3rd April 2023 2:06:34 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Solves an LPN

from svinterface.core.zerod.solver import Solver0Dcpp
from svinterface.core.zerod.lpn import LPN
from svinterface.manager.baseManager import Manager

from pathlib import Path
import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'Solves an LPN')
    
    parser.add_argument('-i', dest = 'config', help = 'path to input lpn')
    parser.add_argument('-sim', type = int, help = 'as_simulation number')
    parser.add_argument('-c', dest = 'csv', action = 'store_true', default = False, help = 'save csv file: Default = False')
    parser.add_argument('-b', dest = 'branch', action = 'store_true', default=False,  help = 'to convert c output to python branch files: Default = False')
    parser.add_argument('--l', dest = 'last_cycle', action = 'store_true', default = False, help = 'only save the last cycle worth of results: Default = False')
    parser.add_argument('--m', dest = 'mean_only', action = 'store_true', default = False, help = 'only save the mean of the results: Default = False')
    parser.add_argument('-v', dest = 'validate', action = 'store_true', default = False, help = 'validate the run with inlet pressure waveform: Default = False')
    
    
    
    args = parser.parse_args()
    
    M = Manager(args.config)
    cur_sim = M['as_simulations'][args.sim]
    
    # get LPN
    lpn_file =  Path(cur_sim['lpn'])
    if not lpn_file.exists():
        raise FileNotFoundError(f'LPN File {str(lpn_file)} not found.')
    lpn = LPN.from_file(str(lpn_file))
    # directory
    rez_dir = Path(cur_sim['dir'])
    # sim
    counter = args.sim
    
    solver = Solver0Dcpp(lpn, last_cycle_only=args.last_cycle, mean_only=args.mean_only, debug = True)
    
    results = solver.run_sim_pipeline(validate = args.validate, save_csv = args.csv, save_branch = args.branch, out_dir = rez_dir)
    
    if args.csv:
        M.register(key = "csv", value = str(rez_dir / "branch_results.csv"), depth = ['as_simulations', counter])    
    if args.branch:
        M.register(key = "npy", value = str(rez_dir / "branch_results.npy"), depth = ['as_simulations', counter])    
    
    M.update()
    
    
    