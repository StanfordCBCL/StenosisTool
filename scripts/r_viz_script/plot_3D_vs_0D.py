
import matplotlib.pyplot as plt
import argparse
import vtk
import numpy as np
from pathlib import Path

from svinterface.core.polydata import Centerlines


def plot_outlets(c_3d: Centerlines, c_1d: Centerlines, save_dir: Path ):
    
    # use valid array
    valid = c_3d.get_pointdata_array("valid")
    outlets = np.where(valid == 1)[0]
    
    # # branch ids
    # br_id = c_3d.get_pointdata_array(c_3d.PointDataFields.BRANCHID)
    # bf_id = c_3d.get_pointdata_array(c_3d.PointDataFields.BIFURCATIONID)

    # # global node id
    # gid = c_3d.get_pointdata_array(c_3d.PointDataFields.NODEID)

    # # outlet points are only connected to one cell (skip inlet point)
    # ids = vtk.vtkIdList()
    # outlets = []
    # for p in range(c_3d.polydata.GetNumberOfPoints()):
    #     c_3d.polydata.GetPointCells(p, ids)
    #     if ids.GetNumberOfIds() == 1 and gid[p] != 0:
    #         assert br_id[p] != -1, 'bifurcation ' + str(bf_id[p]) + ' is connected to an outlet'
    #         outlets.append(gid[p])
    
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
        fig, ax = plt.subplots(1,1)
        ax.plot(results_3d[i]['time'], results_3d[i]['pressure'], label = '3d')
        ax.plot(results_1d[i]['time'], results_1d[i]['pressure'], label = '1d')
        ax.set_ylabel("Pressure (mmHg)")
        ax.set_xlabel("Time (seconds)")
        fig.suptitle(f"Cap {i}: Point {results_3d[i]['point_id']}")
        fig.legend()

        # save
        fig.savefig(str(save_dir / f"cap_{i}.png"))
        plt.close(fig)
    
    
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
        
    
    

if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description = "Plots the differences between 3D and 0D")
    parser.add_argument("-3d", dest = 'three_d', help = '3d results in centerline form (caps only is fine).')
    parser.add_argument("-0d", dest = 'zero_d', help = '-0d results in centerline form (caps only is fine).')
    args = parser.parse_args()
    
    
    c_3d = Centerlines.load_polydata(args.three_d)
    
    c_1d = Centerlines.load_polydata(args.zero_d)
    
    tmp = Path('data/healthy/0080_0001/results/0080_0001/LPN_DIR/0080_0001.sim.0/3D_vs_1D/')
    tmp.mkdir(parents=True, exist_ok = True)
   
    plot_outlets(c_3d, c_1d, tmp)