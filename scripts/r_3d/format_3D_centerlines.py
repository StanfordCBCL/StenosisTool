# File: format_3D_centerlines.py
# File Created: Thursday, 26th January 2023 5:29:38 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Saturday, 11th March 2023 12:21:17 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Converts 3D extracted centerlines to match the 0D extracted form. Assumes extracted centerlines contains exactly 1 time cycle.

from svinterface.core.polydata import Polydata
from svinterface.core.bc import Inflow
from svinterface.manager.baseManager import Manager
import argparse
import re
import numpy as np
from svinterface.utils.misc import d2m

def modify_centerlines(centerlines: Polydata, inflow: Inflow):
    '''maps ts to actual times, as well as change pressures to mmHg
    '''
    
    ts = []
    for arr_name in centerlines.get_pointdata_arraynames():
        x = re.search("pressure_([0-9]*)", arr_name)
        if x:
            ts.append(int(x[1]))
    
    ts = np.array(ts)
    # smooth flow
    inflow.smooth_flow(ts[-1] - ts[0] + 1)
    time = inflow.t[ts - ts[0]]


    
    # fix each one
    
    for f in ['pressure', 'velocity']:
        array_f = []
        if f == 'pressure':
            func = d2m
        else:
            func = lambda x: x
        
        for arr_name in centerlines.get_pointdata_arraynames():
            
            x = re.search(f"{f}_([0-9]+)", arr_name)
            #! assumes in order
            if x:
                array_f.append(func(centerlines.get_pointdata_array(arr_name)))
                cur_ts = int(x[1])
                centerlines.add_pointdata(func(centerlines.get_pointdata_array(arr_name)), arr_name)
                if f == 'velocity':
                    centerlines.rename_pointdata_array(arr_name, f"flow_{inflow.t[cur_ts - ts[0]]:.5f}")
                else:
                    centerlines.rename_pointdata_array(arr_name, f"{f}_{inflow.t[cur_ts - ts[0]]:.5f}")

        array_f = np.array(array_f)
        # compute summary statistics
        print(array_f.shape)
        # avg
        avg = np.trapz(array_f, time, axis = 0) / (time[-1] - time[0])
        centerlines.add_pointdata(avg, 'avg_' + f)
        diff = abs(time - inflow.max_inflow_t)
        # systolic
        sys_tidx = np.argmax(array_f[:, 0])
        #sys_tidx = np.where(diff == min(diff))[0][0]
        centerlines.add_pointdata(array_f[sys_tidx], 'sys_' + f + f'_{time[sys_tidx]:.5f}')
        diff = abs(time - inflow.min_inflow_t)
        # diastolic
        dia_tidx = np.argmin(array_f[:, 0])
        #dia_tidx = np.where(diff == min(diff))[0][0]
        centerlines.add_pointdata(array_f[dia_tidx], 'dia_' + f + f'_{time[dia_tidx]:.5f}')
        
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Converts 3D extracted centerlines to match the 0D extracted form. Assumes extracted centerlines contains exactly 1 time cycle.")
    
    parser.add_argument("-i", dest = 'config', help= 'Config file for project')
    parser.add_argument("-c", dest = 'centerlines', help = "3D extracted centerlines")
    parser.add_argument("-f", dest = "inflow", help = '3D inflow used to compute simulation')
    parser.add_argument("-o", help = "output destination")
    
    args = parser.parse_args()
    
    # manager
    M = Manager(args.config)
    
    
    # centerlines
    c = Polydata.load_polydata(args.centerlines)
    # load inflow
    inflow = Inflow.from_file(args.inflow, inverse = True, smooth = False)
    
    
    modify_centerlines(c, inflow)
    
    c.write_polydata(args.o)
    M.register("3D", args.o, depth = ['workspace'])
    M.update()