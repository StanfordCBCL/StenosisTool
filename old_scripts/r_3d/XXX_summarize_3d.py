# File: summarize_3d.py
# File Created: Sunday, 24th July 2022 1:39:38 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 23rd January 2023 7:42:06 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
#! Description: Turns a 3D simulation with array timesteps into arrays with actual time. Furthermore, summarizes mean, sys and dia, pressures



import numpy as np
import os


import matplotlib.pyplot as plt
from src.misc import create_tool_parser, d2m, get_basename, get_solver_name
from src.polydata import Centerlines
from src.flow import Inflow0D

#########
# Funcs #
#########

def convert_time_steps(centerline3d: Centerlines, tc):
    
    # extract all array names with pressure and velocity
    array_names = centerline3d.get_pointdata_arraynames()
    
    # extract timesteps
    time_steps = []
    for name in array_names:
        if name.startswith('pressure'):
            time_steps.append(int(name.split('_')[-1]))
    
    # get ts per cycle to compute ratio of time step to actual time
    start = min(time_steps) 
    end = max(time_steps)
    ts_per_cycles = end - start
    ts_size = tc/ts_per_cycles
    
    # convert pressures
    for name in array_names:
        if name.startswith('pressure') :
            timestep = int(name.split('_')[-1])
            centerline3d.rename_pointdata(name, 'pressure_' + str(round((timestep - start)*ts_size, 5)) )
        elif name.startswith('velocity'):
            timestep = int(name.split('_')[-1])
            centerline3d.rename_pointdata(name, 'flow_' + str(round((timestep - start)*ts_size, 5)) )
    return
            
def plot_inlet_PAP(centerlines: Centerlines, outfile):
    
    inlet_p = []
    time = []
    array_names = centerlines.get_pointdata_arraynames()
    for name in array_names:
        if name.startswith('pressure'):
            time.append(float(name.split('_')[-1]))
            inlet_p.append(centerlines.get_pointdata(name)[0])
    
    fig, ax = plt.subplots()
    ax.plot(time, inlet_p )
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Inlet Pressure (mmHg)')
    
    mpap = np.trapz(np.array(inlet_p), np.array(time)) / (time[-1] - time[0])
    ax.hlines(y = mpap, xmin = time[0], xmax = time[-1], linewidth=1, color='r', label = 'Avg PAP: ' + str(round(mpap, 2)))
    ax.legend()
    fig.savefig(outfile)
            
def add_summary(centerlines: Centerlines):
    ''' adds max, min, and avg'''
    
    array_names = centerlines.get_pointdata_arraynames()
    
    mpap = {}
    mq = {}
    for name in array_names:
        if name.startswith('pressure'):
            time = float(name.split('_')[-1])
            mpap[time] = centerlines.get_pointdata(name)
            # convert to mmHg
            centerlines.add_pointdata(d2m(mpap[time]), name)
        elif name.startswith('flow'):
            time = float(name.split('_')[-1])
            mq[time] = centerlines.get_pointdata(name)
        
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([d2m(mpap[tidx]) for tidx in t])
    Q = np.array([mq[tidx] for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    mQ = np.trapz(Q, t, axis = 0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'meanPAP')
    centerlines.add_pointdata(mQ, 'meanQ')
    
    
    # find max and min are computed using maxQt/minQt of inlet
    inlet_flows = np.array([mq[tidx][0] for tidx in t])
    t_max = t[np.where(inlet_flows == inlet_flows.max())][0]
    t_min = t[np.where(inlet_flows == inlet_flows.min())][0]

    centerlines.add_pointdata(d2m(mpap[t_max]), 'sysPAP_' + str(t_max))
    centerlines.add_pointdata(d2m(mpap[t_min]), 'diaPAP_' + str(t_min))
    
    centerlines.add_pointdata(mq[t_max], 'sysQ_' + str(t_max))
    centerlines.add_pointdata(mq[t_min], 'diaQ_' + str(t_min))
    
    return centerlines
    


########
# Main #
########

def main(args):
    
    for solver_dir in args.solver_dirs:
        
        
        
        
        join = lambda file: os.path.join(solver_dir, file)
        three_d_dir = join('three_d_dir')
        join_3d = lambda file: os.path.join(three_d_dir, file)
        if not os.path.exists(three_d_dir):
            print(f'{three_d_dir} does not exist. Skipping...')
            continue
        
        # get 3D projected to centerlines vtp file
        solver_file = get_solver_name(solver_dir)
        threed_model_name = get_basename(solver_file)
        three_d_file = join_3d(threed_model_name + '.vtp')
        print(f'Converting 3D file {get_basename(three_d_file)}.')
        if not os.path.exists(three_d_file):
            print(f'{three_d_file} does not exist. Skipping...')
            continue
        
        # get converted filename and check if one already exists
        converted_centerlines = join_3d(threed_model_name + '_converted.vtp')
        if not args.force and os.path.exists(converted_centerlines):
            print(f'A conversion {os.path.basename(converted_centerlines)} already exists. Skipping...')
            continue
            
        
        # load inflow for tc
        inflow_file = join('inflow.flow')
        inflow = Inflow0D.from_file(inflow_file)
        tc = inflow.tc
    
        # load 3D centerlines
        centerlines = Centerlines()
        centerlines.load_centerlines(three_d_file)
        
        
    
        # convert time steps
        print('Converting time steps...', end = '\t', flush = True)
        convert_time_steps(centerline3d=centerlines, tc = tc)
        print('Done')
        print('Adding summary data...', end = '\t', flush = True)
        add_summary(centerlines=centerlines)
        print('Done')
        
        # converted 3D file
        print('Writing centerlines...', end = '\t', flush = True)
        centerlines.write_centerlines(converted_centerlines)
        print('Done')
        
        # plot inlet PAP
        print('Plotting inlet pressures...', end = '\t', flush = True)
        plot_inlet_PAP(centerlines, join_3d('inlet_pressures.png') )
        print('Done')
        



if __name__ == '__main__':
    
    tool = create_tool_parser(desc = 'Converts existing 3D results to summary statistics')
    tool.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'Force overwrite the old file')
    args = tool.parse_args()
    
    main(args)
    