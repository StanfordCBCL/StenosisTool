# File: run_lpn.py
# File Created: Thursday, 3rd November 2022 12:44:51 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 3rd April 2023 1:47:32 pm
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
    parser.add_argument('-sim', type = int, help = 'simulation number')
    parser.add_argument("-n", dest = 'name', default = None, help = 'name of simulation')
    parser.add_argument('-c', dest = 'csv', action = 'store_true', default = False, help = 'save csv file: Default = False')
    parser.add_argument('-b', dest = 'branch', action = 'store_true', default=False,  help = 'to convert c output to python branch files: Default = False')
    parser.add_argument('--l', dest = 'last_cycle', action = 'store_true', default = False, help = 'only save the last cycle worth of results: Default = False')
    parser.add_argument('--m', dest = 'mean_only', action = 'store_true', default = False, help = 'only save the mean of the results: Default = False')
    parser.add_argument('-v', dest = 'validate', action = 'store_true', default = False, help = 'validate the run with inlet pressure waveform: Default = False')
    
    
    
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    # check if simulation exists
    if 'simulations' not in M.yaml or M['simulations'] is None or args.sim not in M['simulations']:
        print("Simulation", args.sim, "was not found. Exitting.")
        exit(1)
        
    sim_data = M['simulations'][args.sim]
    
    # get lpn in simulation
    lpn_file =  Path(sim_data['lpn'])
    if not lpn_file.exists():
        raise FileNotFoundError(f'LPN File {str(lpn_file)} not found.')

    lpn = LPN.from_file(str(lpn_file))
    rez_dir = Path(sim_data['dir'])
    counter = args.sim
    name = args.name if args.name else sim_data['name']
    
    # remove old values & unregister everything
    for file in rez_dir.iterdir():
        file.unlink()
    M.unregister(key=args.sim, depth=['simulations'])
    M.register(key=args.sim,value={}, depth=['simulations'])
    
    
    solver = Solver0Dcpp(lpn, last_cycle_only=args.last_cycle, mean_only=args.mean_only, debug = True)
    
    results = solver.run_sim_pipeline(validate = args.validate, save_csv = args.csv, save_branch = args.branch, out_dir = rez_dir)
    M.register(key = "dir", value = str(rez_dir), depth = ['simulations', counter])
    
    if args.csv:
        M.register(key = "csv", value = str(rez_dir / "branch_results.csv"), depth = ['simulations', counter])    
    if args.branch:
        M.register(key = "npy", value = str(rez_dir / "branch_results.npy"), depth = ['simulations', counter])    
    
    # save a copy of the lpn
    lpn.write_lpn_file(str(rez_dir / lpn_file.name))
    M.register(key="lpn", value = str(rez_dir / lpn_file.name), depth = ['simulations', counter])
    
    M.register(key="name", value = name, depth = ['simulations', counter])
    
    
    M.update()
    
    
    