
from src.solver_results import SolverResults
from src.solver import Solver0D
from src.misc import get_basename, d2m
import sys
import numpy as np
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':

    
    solver_file = sys.argv[1]
    solver = Solver0D()
    solver.read_solver_file(solver_file)
    
    vessel_info = {}
    tree = solver.get_vessel_tree()
    for node in solver.tree_bfs_iterator(tree):
        outlet = False
        if 'boundary_conditions' in node.vessel_info:
            outlet = True
        vessel_info[node.vess_id[0]] = [node.generation, node.side, outlet]
    
    solver_rez_file = os.path.join(os.path.dirname(solver_file), get_basename(solver_file) + '_branch_results.csv')
    solver_rez = SolverResults.load_from_csv(solver_rez_file)
    vessel_names = solver_rez.get_vessel_list()
    for vname in sorted(vessel_names, key = lambda x: int(x[1:])):
            
        vess_df = solver_rez.only_last_cycle(solver_rez.vessel_df(vname), tc = solver.inflow.tc)
        time = np.array(vess_df['time'])
        flow_in = np.array(vess_df['flow_in'])
        flow_out = np.array(vess_df['flow_out'])
        
        dQ = (flow_in - flow_out)/flow_in
        pressure_in = np.array(vess_df['pressure_in'])
        pressure_out = np.array(vess_df['pressure_out'])
        dp = pressure_in - pressure_out
        mpd = np.trapz(dp, time)/solver.inflow.tc
        mpoutd = np.trapz(pressure_out, time)/solver.inflow.tc
        
        
        '''
        fig, ax = plt.subplots(3, 1, sharex = True)
        ax[0].plot(time, d2m(pressure_in - pressure_out))
        ax[0].set_ylabel('pressure_diff (mmHg)')
        ax[1].plot(time, d2m( pressure_out))
        ax[1].set_ylabel('pressure_out (mmHg)')
        ax[2].plot(time, d2m( pressure_in))
        ax[2].set_ylabel('pressure_in (mmHg)')
        ax[2].set_xlabel('time')
        fig.savefig(f'scratch/{vname}_pressure.png')
        plt.close(fig)
        fig, ax = plt.subplots(2, 1,sharex = True)
        ax[0].plot(time, flow_in)
        ax[0].set_ylabel('flow_in (mL/s)')
        ax[1].plot(time, flow_out)
        ax[1].set_xlabel('time')
        ax[1].set_ylabel('flow_out (mL/s)')
        fig.savefig(f'scratch/{vname}_flows.png')
        plt.close(fig)
        '''
        mdQ = np.trapz(dQ,  time)/solver.inflow.tc
        print(vname + ': meanPdrop =', round(d2m(mpd),2) , 'meanPout =', round(d2m(mpoutd),2) , 'mQdrop =', mdQ , 'Generation = ', vessel_info[int(vname[1:])][0], 'Side =', vessel_info[int(vname[1:])][1], 'Outlet =', vessel_info[int(vname[1:])][2])
    
        
    