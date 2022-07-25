
# TODO: Summarize the 3D model results into usable summary arrays usable for error calculation 
from audioop import add
import numpy as np
import shutil
import os


import matplotlib.pyplot as plt

from src.data_org import DataPath, JunctionStenosisResults
from src.misc import create_parser, d2m
from src.centerlines import Centerlines

#########
# Funcs #
#########

def convert_time_steps(centerline3d: Centerlines, params3d):
    
    # extract all array names with pressure and velocity
    # ! I DO NOT KNOW IF VELOCITY = FLOW so start with pressures only
    array_names = centerline3d.get_pointdata_arraynames()
    
    # get number of timesteps per cycle, and extract last cycle
    ts_per_cycles = int(params3d['sim_steps_per_cycle'])
    tc = float(params3d['sim_period'])
    ts_start = int(params3d['sim_vtu_start'])
    
    # compute ratio of timestep to actual time
    ts_size = tc/ts_per_cycles
    
    # extract timesteps and compute last cycle
    time_steps = []
    for name in array_names:
        if name.startswith('pressure'):
            time_steps.append(int(name.split('_')[-1]))
    
    last_cycle_start = max(time_steps) - ts_per_cycles - ts_start
    last_cycle_end = max(time_steps) - ts_start

    # convert pressures
    pressures = {}
    for name in array_names:
        if name.startswith('pressure'):
            timestep = int(name.split('_')[-1])
            if timestep < last_cycle_start or timestep > last_cycle_end:
                centerline3d.remove_pointdata(name)
            else:
                centerline3d.rename_pointdata(name, 'pressure_' + str(round((timestep - last_cycle_start)*ts_size, 5)) )
        #! Temporary
        if name.startswith('velocity'):
            centerline3d.remove_pointdata(name)

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
            
def add_summary_PAP(centerlines: Centerlines):
    ''' adds max, min, and avg'''
    
    array_names = centerlines.get_pointdata_arraynames()
    mpap = {}

    for name in array_names:
        if name.startswith('pressure'):
            time = float(name.split('_')[-1])
            mpap[time] = centerlines.get_pointdata(name)
            # convert to mmHg
            centerlines.add_pointdata(d2m(mpap[time]), name)
            

            

            
    
    t = np.array(list(sorted(mpap.keys())))
    P = np.array([d2m(mpap[tidx]) for tidx in t])
    mPAP = np.trapz(P, t, axis = 0) / (t[-1] - t[0])
    centerlines.add_pointdata(mPAP,'mPAP')
    
    
    # Rough Approximation of max and min PAP
    cur_min = 0
    cur_max = 0
    for time in t[1:]:
        if mpap[time].max() < mpap[cur_min].max():
            cur_min = time
        if mpap[time].max() > mpap[cur_max].max():
            cur_max = time

    centerlines.add_pointdata(d2m(mpap[cur_max]), 'maxPAP_' + str(cur_max))
    centerlines.add_pointdata(d2m(mpap[cur_min]), 'minPAP_' + str(cur_min))
    
    return centerlines
    


########
# Main #
########

def tool_main(args):
    raise NotImplementedError

def dev_main(args):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        model = org.find_model(model_name)
        
        model_results = JunctionStenosisResults(args.root, model)

        if model_name not in model_results.available_3d_models:
            print('3D Model is not available. Skipping.')
            continue
        
        # copy over centerlines file
        new_centerlines_file = os.path.join(model_results.model_three_d_dir, model_name + '.vtp')
        shutil.copy(os.path.join(model_results.three_d_reruns, model_name + '.vtp'), model_results.model_three_d_dir)
        
        # parameters
        model_3d_params = np.load(model_results.three_d_parameters, allow_pickle=True).item()[model_name]['params']
        
        # load centerlines
        centerlines = Centerlines()
        centerlines.load_centerlines(new_centerlines_file)
        
        # convert time steps
        
        convert_time_steps(centerline3d=centerlines, params3d=model_3d_params)
        add_summary_PAP(centerlines=centerlines)
        
        centerlines.write_centerlines(new_centerlines_file)
        
        plot_inlet_PAP(centerlines, os.path.join(model_results.model_three_d_dir,'inlet_pressures.png') )
        



if __name__ == '__main__':
    
    parser, dev, tool = create_parser(desc = 'Converts existing 3D results to summary statistics')
    
    
    args = parser.parse_args()
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)