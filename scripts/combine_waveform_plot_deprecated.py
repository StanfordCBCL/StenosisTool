# File: combine_waveform_plot_deprecated.py
# File Created: Wednesday, 27th July 2022 9:37:19 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 5th August 2022 1:37:29 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Can combine waveforms given a multitude of csv result files. Deprecated due to difficulty of use, but technically works


from src.solver import Solver0D
from src.solver_results import SolverResults
from src.misc import d2m, get_basename
import matplotlib.pyplot as plt
import numpy as np
import argparse
def combine_solver_results(solver: Solver0D, solver_results: dict, out_file):
    # retrieve time of cardiac cycle and number of points per cycle
    inflow_tc = solver.inflow.tc
    num_pts = int(solver.simulation_params['number_of_time_pts_per_cardiac_cycle'])
        
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(1,1 ,figsize=(15, 10))
    ax.set_title(f"Inlet Pressure", fontdict={'fontsize': 24})
    ax.set_xlabel('time (s)', fontdict={'fontsize': 20})
    ax.set_ylabel('pressure (mmHg)', fontdict={'fontsize': 20})
    # plot last 3 cycles
    last_three_cycles = -3 *  num_pts
    
    for name, sim_results in solver_results.items():
        v0 = sim_results.vessel_df('V0')
        inlet_pressure = d2m(np.array(v0['pressure_in'][last_three_cycles:]))
        time_ltc = np.array(v0['time'][last_three_cycles:])
        # plot pressure curves
        ax.plot(time_ltc, inlet_pressure, label = name)
        # get last cardiac cycle
        v0_lc = SolverResults.only_last_cycle(v0, tc = inflow_tc)
        
    ax.tick_params(axis="x", labelsize=16) 
    ax.tick_params(axis = 'y', labelsize=16)
    ax.legend(fontsize = 24, loc = 'upper left', framealpha = .5)

    fig.savefig(out_file)
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', dest = 'solver_file', help = 'base solver file to retireve tc from')
    parser.add_argument('-r', dest = 'result_files',action = 'append', default = [], nargs = 2, help = 'Results/name files to combine')
    parser.add_argument('-o', dest = 'outfile', help = 'outfile')
    
    args = parser.parse_args()
    solver_results = {}
    for result_file, name in args.result_files:
        solver_results[name] = SolverResults.load_from_csv(result_file)
    
    solver = Solver0D()
    print(args.solver_file)
    solver.read_solver_file(args.solver_file)
    combine_solver_results(solver, solver_results, out_file = args.outfile)