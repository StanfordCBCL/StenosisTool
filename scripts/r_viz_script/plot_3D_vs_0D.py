
import matplotlib.pyplot as plt
import argparse
import vtk
import numpy as np
from pathlib import Path


from svinterface.core.polydata import Centerlines
from svinterface.core.zerod.solver import SolverResults
from svinterface.manager.baseManager import Manager
from svinterface.utils.io import write_json


def plot_outlets(c_3d: Centerlines, c_1d: Centerlines, save_dir: Path ):
    
    # use valid array
    valid = c_3d.get_pointdata_array("Caps_0D")
    outlets = np.where(valid != -1)[0]
    outlets = np.array([0] + list(outlets))
    
    results_3d = {}
    # iterate through each outlet
    for oidx, point_id in enumerate(outlets):
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
    for oidx, point_id in enumerate(outlets):
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
    # save each outlet
    for i in range(len(outlets)):
        if i == 0:
            print(max(results_3d[i]['pressure']))
        fig, ax = plt.subplots(1,1)
        ax.plot(results_3d[i]['time'], results_3d[i]['pressure'], label = '3d')
        ax.plot(results_1d[i]['time'], results_1d[i]['pressure'], label = '0d')
        ax.set_ylabel("Pressure (mmHg)")
        ax.set_xlabel("Time (seconds)")
        fig.suptitle(f"Cap {i}: Point {results_3d[i]['point_id']}")
        fig.legend()

        # save
        fig.savefig(str(save_dir / f"cap_{i}.png"))
        plt.close(fig)
    
    
    errors = {}
    # means
    fig, ax = plt.subplots(1, 1)
    oned_means = []
    threed_means = []
    for i in range(len(outlets)):
        oned_means.append(np.trapz(results_1d[i]['pressure'], results_1d[i]['time']) / (results_1d[i]['time'][-1] - results_1d[i]['time'][0]))
        threed_means.append(np.trapz(results_3d[i]['pressure'], results_3d[i]['time']) / (results_3d[i]['time'][-1] - results_3d[i]['time'][0]))
    ax.scatter(range(len(outlets)), threed_means, label = '3d')
    ax.scatter(range(len(outlets)), oned_means, label = '0d')
    fig.legend()
    fig.suptitle("Mean of outlets")
    ax.set_ylabel("Mean Pressure (mmHg)")
    ax.set_xlabel("Outlets")
    
    fig.savefig(str(save_dir / "outlet_means.png"))
    plt.close(fig)
    
    errors_means = (np.array(threed_means) - np.array(oned_means))
    errors['Mean'] = {}
    errors['Mean']['RootMSE'] = np.sqrt((errors_means ** 2).mean()).item()
    errors['Mean']['MeanAE'] = np.abs(errors_means).mean().item()
    errors['Mean']['MaxAE'] = np.abs(errors_means).max().item()
    errors['Mean']['MinAE'] = np.abs(errors_means).min().item()
    errors['Mean']['MeanRelativeError'] = np.abs(errors_means / np.array(threed_means)).mean().item()
    
    
    # systolic
    fig, ax = plt.subplots(1, 1)
    oned_maxs = []
    threed_maxs = []
    for i in range(len(outlets)):
        oned_maxs.append(np.array(results_1d[i]['pressure']).max())
        threed_maxs.append(np.array(results_3d[i]['pressure']).max())
    ax.scatter(range(len(outlets)), threed_maxs, label = '3d')
    ax.scatter(range(len(outlets)), oned_maxs, label = '0d')
    fig.legend()
    fig.suptitle("Systolic pressures of outlets")
    ax.set_ylabel("Systolic Pressure (mmHg)")
    ax.set_xlabel("Outlets")
    
    fig.savefig(str(save_dir / "outlet_maxs.png"))
    plt.close(fig)
    
    errors_maxs = (np.array(threed_maxs) - np.array(oned_maxs))
    errors['Max'] = {}
    errors['Max']['RootMSE'] = np.sqrt((errors_maxs ** 2).mean()).item()
    errors['Max']['MeanAE'] = np.abs(errors_maxs).mean().item()
    errors['Max']['MaxAE'] = np.abs(errors_maxs).max().item()
    errors['Max']['MinAE'] = np.abs(errors_maxs).min().item()
    errors['Max']['MeanRelativeError'] = np.abs(errors_maxs / np.array(threed_maxs)).mean().item()
    
    # diastolic
    fig, ax = plt.subplots(1, 1)
    oned_mins = []
    threed_mins = []
    for i in range(len(outlets)):
        oned_mins.append(np.array(results_1d[i]['pressure']).min())
        threed_mins.append(np.array(results_3d[i]['pressure']).min())
    ax.scatter(range(len(outlets)), threed_mins, label = '3d')
    ax.scatter(range(len(outlets)), oned_mins, label = '0d')
    fig.legend()
    fig.suptitle("Diastolic pressures of outlets")
    ax.set_ylabel("Diastolic Pressure (mmHg)")
    ax.set_xlabel("Outlets")
    
    fig.savefig(str(save_dir / "outlet_mins.png"))
    plt.close(fig)
    
    errors_mins = (np.array(threed_mins) - np.array(oned_mins))
    errors['Min'] = {}
    errors['Min']['RootMSE'] = np.sqrt((errors_mins ** 2).mean()).item()
    errors['Min']['MeanAE'] = np.abs(errors_mins).mean().item()
    errors['Min']['MaxAE'] = np.abs(errors_mins).max().item()
    errors['Min']['MinAE'] = np.abs(errors_mins).min().item()
    errors['Min']['MeanRelativeError'] = np.abs(errors_mins / np.array(threed_mins)).mean().item()
    
    # errors for MPA
    errors['MPA'] = {'Sys': {}, 'Mean': {}, 'Dia': {}}
    errors['MPA']['Sys']['AbsoluteError'] = abs(threed_maxs[0] - oned_maxs[0])
    errors['MPA']['Sys']['RelativeError'] = (threed_maxs[0] - oned_maxs[0]) / threed_maxs[0]
    errors['MPA']['Mean']['AbsoluteError'] = abs(threed_means[0] - oned_means[0])
    errors['MPA']['Mean']['RelativeError'] = (threed_means[0] - oned_means[0]) / threed_means[0]
    errors['MPA']['Dia']['AbsoluteError'] = abs(threed_mins[0] - oned_mins[0])
    errors['MPA']['Dia']['RelativeError'] = (threed_mins[0] - oned_mins[0]) / threed_mins[0]
    
    write_json(save_dir / "errors.json", errors)
    
    # # plot distribution of differences
    oned_means = np.array(oned_means)
    threed_means = np.array(threed_means)
    
    fig, ax = plt.subplots(1, 1)
    ax.hist((threed_means - oned_means )/ threed_means, bins = 'auto')
    fig.savefig(str(save_dir / "outlet_diffs.png"))
    plt.close(fig)
    
    
    

if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description = "Plots the differences between 3D and 0D")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file (must contain a simulation)')
    parser.add_argument('-mode', default = None, help = 'Mode: None, AS, R')
    parser.add_argument("-sim", type = int, dest = 'sim', help = '0d simulation results desired to compare')
    parser.add_argument("--caps", action = "store_true", default = False, help = "caps only")
    parser.add_argument("--junctions", action = "store_true", default = False, help = "junctions only")
    
    
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    sims = 'simulations'
    if args.mode == 'AS':
        sims = 'as_simulations'
    elif args.mode == 'R':
        sims = 'r_simulations'
    # else:
    #     raise ValueError("-mode must be AS or R or not set")
        
    
    three_d = M['workspace']['3D']
    try:
        zero_d_sim = Path(M[sims][args.sim]['dir'])
    except IndexError:
        print(f"Simulation {args.sim} not found")
        exit(1)
    
    c_3d = Centerlines.load_polydata(three_d)
    
    csv_0d = zero_d_sim / "branch_results.csv"
    if not csv_0d.exists():
        print(f"Simulation {args.sim} does not contain a branch_results.csv file.")
        exit(1)
    r_0d = SolverResults.from_csv(str(csv_0d))
    
    c_0d = Centerlines.load_centerlines(M[sims][args.sim]['centerlines'])
    
    comp_folder = zero_d_sim / "3D_vs_0D"
    comp_folder.mkdir(exist_ok=True, parents=True)
   
    if args.caps:
        cap_folder = comp_folder / "caps"
        cap_folder.mkdir(exist_ok=True, parents=True)
        plot_outlets(c_3d, c_0d, cap_folder)
        
    # if args.junctions:
    #     j_folder = comp_folder / "junctions"
    #     j_folder.mkdir(exist_ok=True, parents=True)
    #     plot_junctions(c_3d, r_0d, j_folder)