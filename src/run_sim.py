
from .solver import Solver0D
from .solver_results import SolverResults
from .misc import *
from svzerodsolver.runnercpp import run_from_config
import matplotlib.pyplot as plt
import numpy as np
import os


def get_branch_results_file(solver_file, cpp = False):
    ''' Finds the branch result file name given the solver file. (cpp or python)
    '''
    if cpp:
        return os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_branch_results.csv'
    else:
        return os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_branch_results.npy'
    
def get_waveform_file(solver_file):
    ''' finds waveform file name given the solver file
    '''
    return os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_inlet_pressures.png'
    

def validate_rez(solver: Solver0D, sim_results: SolverResults, out_file):
    ''' saves inlet waveforms as solver_file_inlet_pressure.png'''
    # retrieve time of cardiac cycle and number of points per cycle
    inflow_tc = solver.inflow.tc
    num_pts = int(solver.simulation_params['number_of_time_pts_per_cardiac_cycle'])
    
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(1,1 ,figsize=(15, 10))

    # plot last 3 cycles
    last_three_cycles = -3 *  num_pts
    v0 = sim_results.vessel_df('V0')
    inlet_pressure = d2m(np.array(v0['pressure_in'][last_three_cycles:]))
    time_ltc = np.array(v0['time'][last_three_cycles:])
    # plot pressure curves
    ax.plot(time_ltc, inlet_pressure)
    ax.set_title(f"Inlet Pressure", fontdict={'fontsize': 24})
    ax.set_xlabel('time (s)', fontdict={'fontsize': 20})
    ax.set_ylabel('pressure (mmHg)', fontdict={'fontsize': 20})
    
    # get last cardiac cycle
    v0_lc = SolverResults.only_last_cycle(v0, tc = inflow_tc)

    # compute mean, max, min PAP
    stable_inlet_pressure = d2m(np.array(v0_lc['pressure_in']))
    time_lc = np.array(v0_lc['time'])
    mpap_sim = np.trapz(stable_inlet_pressure, time_lc) / inflow_tc
    
    x_max = time_lc[np.where(stable_inlet_pressure == stable_inlet_pressure.max())] + time_ltc[0] + 2 * inflow_tc

    x_min = time_lc[np.where(stable_inlet_pressure == stable_inlet_pressure.min())] + time_ltc[0] + 2 * inflow_tc
    
    # plot mean, max, min
    ax.plot(x_max, 
            stable_inlet_pressure.max(), 
            'r^', 
            label = 'Max PAP: ' + str(round(stable_inlet_pressure.max(),2)))
    ax.plot(x_min, 
            stable_inlet_pressure.min(), 
            'g^', 
            label = 'Min PAP: ' + str(round(stable_inlet_pressure.min(), 2)))
    ax.hlines(y = mpap_sim, xmin = time_ltc[0], xmax = time_ltc[-1], linewidth=1, color='b', label = 'Avg PAP: ' + str(round(mpap_sim, 2)) )
    
    ax.tick_params(axis="x", labelsize=16) 
    ax.tick_params(axis = 'y', labelsize=16)
    ax.legend(fontsize = 24, loc = 'upper left', framealpha = .5)

    fig.savefig(out_file)
    
    
def run_sim(solver: Solver0D, use_steady_soltns = True, save_branch_results = False, save_csv = False, save_last = False,  debug = False):
    ''' Takes in solver file and runs a simulation'''
    # Run a c simulation to test
    
    # disable steady solutions if requested
    if not use_steady_soltns:
        solver.simulation_params['steady_initial'] = use_steady_soltns
    
    # run cpp
    if debug:
        print('Running Solver...', end = '\t', flush = True)
    results = run_from_config(solver.solver_data)
    rez = SolverResults(results)
    if debug:
        print('Done')
        
    if save_last:
        rez_out = SolverResults(results, last_cycle=True, tc = solver.inflow.tc)
        rez_out.result_df = SolverResults.only_last_cycle(rez_out.result_df, tc = solver.inflow.tc)
    else:
        rez_out = rez
        
    if save_csv:
        if debug:
            print('Saving CSV...', end = '\t', flush = True)
        rez_out.save_csv(get_branch_results_file(solver.solver_file, cpp = True))
        if debug:
            print('Done')
    if save_branch_results:
        if debug:
            print('Converting to branch results...', end = '\t', flush = True)
        rez_out.convert_to_python(solver, get_branch_results_file(solver.solver_file, cpp = False))
        if debug:
            print('Done')
    return rez
