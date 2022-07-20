from .file_io import *
from .misc import *
import svzerodsolver as zerod
import matplotlib.pyplot as plt
import numpy as np
import os

def get_result_file(solver_file):
    ''' Finds the result file name given the solver file.'''
    return os.path.join(os.path.dirname(solver_file),os.path.splitext(os.path.basename(solver_file))[0]) + '_all_results.npy'


def validate_rez(solver_file, waveform_name):
    ''' saves inlet waveforms as solver_file_inlet_pressure.png'''
    
    solver = Solver0D()
    solver.read_solver_file(solver_file)
    for bc in solver.bc:
        if bc['bc_name'] == 'INFLOW':
            break
    # retrieve time of cardiac cycle and number of points per cycle
    inflow_tc = bc['bc_values']['t'][-1]
    num_pts = int(solver.simulation_params['number_of_time_pts_per_cardiac_cycle'])
        
    sim_results = SolverResults(get_result_file(solver_file))
    
    ## save flow and pressure graphs (last 3 cycles)
    fig, ax = plt.subplots(1,1 ,figsize=(15, 10))

    # plot last 3 cycles
    last_three_cycles = -3 *  num_pts
    inlet_pressure = d2m(sim_results.pressures['P_BC0_inlet_V0'][last_three_cycles:])
    # plot pressure curves
    ax.plot(sim_results.time[last_three_cycles:], inlet_pressure)
    ax.set_title(f"Inlet Pressure")
    ax.set_xlabel('time (s)')
    ax.set_ylabel('pressure (mmHg)')
    
    # get last cardiac cycle
    last_cycle = -1 * num_pts
    time = sim_results.time[last_cycle:]
    time = time - time[0] + sim_results.time[0]

    # compute mean, max, min PAP
    stable_inlet_pressure = d2m(sim_results.pressures['P_BC0_inlet_V0'][last_cycle:])
    mpap_sim = np.trapz(stable_inlet_pressure, time) / inflow_tc
    
    x_max = sim_results.time[last_cycle:][np.where(stable_inlet_pressure == stable_inlet_pressure.max())]

    
    x_min = sim_results.time[last_cycle:][np.where(stable_inlet_pressure == stable_inlet_pressure.min())]
    
    # plot mean, max, min
    ax.plot(x_max, 
            stable_inlet_pressure.max(), 
            'r^', 
            label = 'Max PAP: ' + str(stable_inlet_pressure.max()))
    ax.plot(x_min, 
            stable_inlet_pressure.min(), 
            'g^', 
            label = 'Min PAP: ' + str(stable_inlet_pressure.min()))
    ax.hlines(y = mpap_sim, xmin = sim_results.time[last_three_cycles], xmax = sim_results.time[-1], linewidth=1, color='b', label = 'Avg PAP: ' + str(mpap_sim) )
    
    ax.legend()

    waveform_name = os.path.splitext(os.path.basename(solver_file))[0] + '_inlet_pressures.png'
    fig.savefig(os.path.join(os.path.dirname(solver_file), waveform_name))
    
    
def run_sim(solver_file, save_branches = True, block_print = True, use_steady_soltns = True, last_cycle = False):
    # Run a simulation to test
    
    if block_print:
        blockPrint() # silence
    zerod.solver.set_up_and_run_0d_simulation(zero_d_solver_input_file_path=solver_file,
                                    last_cycle=last_cycle,
                                    save_results_all=True,
                                    save_results_branch=save_branches,
                                    use_steady_soltns_as_ics = use_steady_soltns)
    if block_print:
        enablePrint()

