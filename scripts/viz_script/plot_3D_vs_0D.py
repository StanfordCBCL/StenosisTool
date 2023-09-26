
import matplotlib.pyplot as plt
import argparse
import vtk
import numpy as np
from pathlib import Path
import shutil

from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.solver import SolverResults
from svinterface.manager.baseManager import Manager
from svinterface.utils.io import write_json

from svinterface.plotting.params import set_params


def plot_valid(c_3d: Centerlines, c_1d: Centerlines, save_dir: Path, points ):
    
    # use valid array
    caps = c_3d.get_pointdata_array("Caps_0D")
    juncs = c_3d.get_pointdata_array("Junctions_0D") 
    vess = c_3d.get_pointdata_array("Vessels_0D") 
    #! pull out which outlet it actually is
    valid = np.array(sorted(list(set([0] + list(np.where(caps != -1)[0]) + list(np.where(juncs != -1)[0]) + list(np.where(vess != -1)[0])))))
    
    results_3d = {}
    # iterate through each valid point
    for oidx, point_id in enumerate(valid):
        time = []
        pressure = []
        flow = []
        #! use the fact that they should be in order already
        for arr_name in c_3d.get_pointdata_arraynames():
            if arr_name.startswith("pressure_"):
                time.append(float(arr_name.split('_')[1]))
                pressure.append(c_3d.polydata.GetPointData().GetArray(arr_name).GetValue(point_id))
        
        results_3d[oidx] = {'time': time,
                            'pressure': pressure,
                            'flow': flow,
                            'point_id': point_id}
        
    results_1d = {}
    # iterate through each outlet
    for oidx, point_id in enumerate(valid):
        time = []
        pressure = []
        flow = []
        #! use the fact that they should be in order already
        for arr_name in c_1d.get_pointdata_arraynames():
            if arr_name.startswith("pressure_"):
                time.append(float(arr_name.split('_')[1]))
                pressure.append(c_1d.polydata.GetPointData().GetArray(arr_name).GetValue(point_id))
        
        results_1d[oidx] = {'time': time,
                            'pressure': pressure,
                            'flow': flow,
                            'point_id': point_id}
    
    # make dir if it doesn't exist
    save_dir.mkdir(parents = True, exist_ok=True)
    
    ## Summary Statistics
    
    errors = {}
    ## means
    fig, ax = plt.subplots(1, 3, figsize=(26, 8))
    fig.suptitle("Summary Comparison of 0D and 3D Values at Relevant Points")
    zerod_means = []
    threed_means = []
    for i in range(len(valid)):
        zerod_means.append(np.trapz(results_1d[i]['pressure'], results_1d[i]['time']) / (results_1d[i]['time'][-1] - results_1d[i]['time'][0]))
        threed_means.append(np.trapz(results_3d[i]['pressure'], results_3d[i]['time']) / (results_3d[i]['time'][-1] - results_3d[i]['time'][0]))
    ax[1].scatter(range(len(valid)), threed_means, label = '3d')
    ax[1].scatter(range(len(valid)), zerod_means, label = '0d')
    ax[1].legend()
    ax[1].set_title("Mean")
    ax[1].set_ylabel("Pressure (mmHg)")
    ax[1].set_xlabel("Points")
    
    errors_means = (np.array(threed_means) - np.array(zerod_means))
    errors['Mean'] = {}
    errors['Mean']['RootMSE'] = np.sqrt((errors_means ** 2).mean()).item()
    errors['Mean']['MeanAE'] = np.abs(errors_means).mean().item()
    errors['Mean']['MaxAE'] = np.abs(errors_means).max().item()
    errors['Mean']['MinAE'] = np.abs(errors_means).min().item()
    errors['Mean']['MeanRelativeError'] = np.abs(errors_means / np.array(threed_means)).mean().item()
    errors['Mean']['MaxRelativeError'] = np.abs(errors_means / np.array(threed_means)).max().item()
    
    ## systolic
    zerod_maxs = []
    threed_maxs = []
    for i in range(len(valid)):
        zerod_maxs.append(np.array(results_1d[i]['pressure']).max())
        threed_maxs.append(np.array(results_3d[i]['pressure']).max())
    ax[0].scatter(range(len(valid)), threed_maxs, label = '3d')
    ax[0].scatter(range(len(valid)), zerod_maxs, label = '0d')
    ax[0].legend()
    ax[0].set_title("Systolic")
    ax[0].set_ylabel("Pressure (mmHg)")
    ax[0].set_xlabel("Points")
    
    errors_maxs = (np.array(threed_maxs) - np.array(zerod_maxs))
    errors['Max'] = {}
    errors['Max']['RootMSE'] = np.sqrt((errors_maxs ** 2).mean()).item()
    errors['Max']['MeanAE'] = np.abs(errors_maxs).mean().item()
    errors['Max']['MaxAE'] = np.abs(errors_maxs).max().item()
    errors['Max']['MinAE'] = np.abs(errors_maxs).min().item()
    errors['Max']['MeanRelativeError'] = np.abs(errors_maxs / np.array(threed_maxs)).mean().item()
    errors['Max']['MaxRelativeError'] = np.abs(errors_maxs / np.array(threed_maxs)).max().item()
    
    # diastolic
    zerod_mins = []
    threed_mins = []
    for i in range(len(valid)):
        zerod_mins.append(np.array(results_1d[i]['pressure']).min())
        threed_mins.append(np.array(results_3d[i]['pressure']).min())
    ax[2].scatter(range(len(valid)), threed_mins, label = '3d')
    ax[2].scatter(range(len(valid)), zerod_mins, label = '0d')
    ax[2].legend()
    ax[2].set_title("Diastolic")
    ax[2].set_ylabel("Pressure (mmHg)")
    ax[2].set_xlabel("Points")
    
    fig.savefig(str(save_dir / "summary.png"))
    plt.close(fig)
    
    errors_mins = (np.array(threed_mins) - np.array(zerod_mins))
    errors['Min'] = {}
    errors['Min']['RootMSE'] = np.sqrt((errors_mins ** 2).mean()).item()
    errors['Min']['MeanAE'] = np.abs(errors_mins).mean().item()
    errors['Min']['MaxAE'] = np.abs(errors_mins).max().item()
    errors['Min']['MinAE'] = np.abs(errors_mins).min().item()
    errors['Min']['MeanRelativeError'] = np.abs(errors_mins / np.array(threed_mins)).mean().item()
    errors['Min']['MaxRelativeError'] = np.abs(errors_mins / np.array(threed_mins)).max().item()
    
    # errors for MPA
    errors['MPA'] = {'Sys': {}, 'Mean': {}, 'Dia': {}}
    errors['MPA']['Sys']['AbsoluteError'] = abs(threed_maxs[0] - zerod_maxs[0])
    errors['MPA']['Sys']['RelativeError'] = (threed_maxs[0] - zerod_maxs[0]) / threed_maxs[0]
    errors['MPA']['Mean']['AbsoluteError'] = abs(threed_means[0] - zerod_means[0])
    errors['MPA']['Mean']['RelativeError'] = (threed_means[0] - zerod_means[0]) / threed_means[0]
    errors['MPA']['Dia']['AbsoluteError'] = abs(threed_mins[0] - zerod_mins[0])
    errors['MPA']['Dia']['RelativeError'] = (threed_mins[0] - zerod_mins[0]) / threed_mins[0]
    
    write_json(save_dir / "errors.json", errors)
    
    ## plot 3D on x axis, and 0D on y axis
    fig, ax = plt.subplots(1, 3, figsize=(26, 8))
    fig.suptitle("3D vs 0D Pressures at Relevant Points.")
    ax[0].plot([min(threed_maxs), max(threed_maxs)], [min(threed_maxs), max(threed_maxs)])
    ax[0].scatter(threed_maxs, zerod_maxs, c="r")
    ax[0].set_title("Systolic")
    ax[1].plot([min(threed_means), max(threed_means)], [min(threed_means), max(threed_means)])
    ax[1].scatter(threed_means, zerod_means, c="r")
    ax[1].set_title("Mean")
    ax[2].plot([min(threed_mins), max(threed_mins)], [min(threed_mins), max(threed_mins)])
    ax[2].scatter(threed_mins, zerod_mins, c="r")
    ax[2].set_title("Diastolic")
    for i in range(3):
        ax[i].set_ylabel("0D Pressure (mmHg)")
        ax[i].set_xlabel("3D Pressure (mmHg)")
    
    
    fig.savefig(str(save_dir / "summary.2.png"))
    plt.close(fig)
    
    if not points:
        return
    
    ## individual outlets
    individual_dir = save_dir / "points"
    individual_dir.mkdir(exist_ok=True)
    # save each outlet
    for i in range(len(valid)):
        if i == 0:
            print(max(results_3d[i]['pressure']))
        fig, ax = plt.subplots(1,1)
        ax.plot(results_3d[i]['time'], results_3d[i]['pressure'], label = '3d')
        ax.plot(results_1d[i]['time'], results_1d[i]['pressure'], label = '0d')
        ax.set_ylabel("Pressure (mmHg)")
        ax.set_xlabel("Time (seconds)")
        fig.suptitle(f"Point {results_3d[i]['point_id']}")
        fig.legend()

        # save
        fig.savefig(str(individual_dir / f"gid_{i}.png"))
        plt.close(fig)
    
    

if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description = "Plots the differences between 3D and 0D")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file (must contain a simulation)')
    parser.add_argument('-mode', default = None, help = 'Mode: None, P')
    parser.add_argument('-3D', dest = 'threed', default = None, required = False, help = '3D formatted centerlines')
    parser.add_argument("-sim", dest = 'sim', help = '0d simulation results desired to compare')
    parser.add_argument('-points', action = 'store_true', default=False, help = 'whether to plot individual points')
    
    
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    if args.mode is None:
        sim = M['simulations'][int(args.sim)]
    elif args.mode == 'P':
        sim = M['parameterization']['corrections'][args.sim]
    
    three_d = M['workspace']['3D'] if (args.threed is None) else args.threed
    try:
        zero_d_sim = Path(sim['dir'])
    except IndexError:
        print(f"Simulation {args.sim} not found")
        exit(1)
    
    c_3d = Centerlines.load_polydata(three_d)
    
    c_0d = Centerlines.load_centerlines(sim['centerlines'])
    
    comp_folder = zero_d_sim / "3D_vs_0D"
    if comp_folder.exists():
        shutil.rmtree(str(comp_folder))
    comp_folder.mkdir(exist_ok=True, parents=True)
    
    # plotting params
    set_params()
    
    plot_valid(c_3d, c_0d, comp_folder, points = args.points)