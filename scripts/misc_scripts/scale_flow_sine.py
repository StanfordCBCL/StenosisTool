# File: scale_flow_sine.py
# File Created: Monday, 22nd May 2023 12:02:17 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 31st May 2023 2:58:18 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Scales inflow to be appropriate using 2 sinusoids. Does not shift systolic

from svinterface.core.bc import Inflow
from svinterface.plotting.plot_flow import plot_flow

import argparse
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Add sine to plot")
    
    parser.add_argument("-i", help = 'inflow waveform')
    parser.add_argument("-o", help = "inflow waveform scaled out file")
    parser.add_argument("--co", dest = "co", type=float, help = "cardiac outflow")
    
    args = parser.parse_args()

    inflow = Inflow.from_file(args.i, smooth=False)
    # plot_flow(inflow)
    
    
    # hard coded points for AS1 flow.
    startpoint = 0
    midpoint = inflow.max_inflow_t
    # endpoint = .5801
    endpoint = inflow.min_inflow_t
    
    
    # compute diff
    print(args.co * 1000/60, inflow.mean_inflow)
    co_diff = args.co * 1000/60 - inflow.mean_inflow
    
    # #####
    # # 4th order
    # #####
    # t_idx = np.where(inflow.t <= endpoint)[0]
    # t = inflow.t[t_idx]
    
    # y_init = (t - midpoint)**2 * (t - endpoint) * t
    # y = (co_diff / np.trapz(y_init, t)) * y_init
    
    # inflow.Q[t_idx] += y
    
    
    #####
    # disjoint Sinusoid
    #####
    
    # split ratio by xrange
    co_diff1 = co_diff * (midpoint - startpoint)/(endpoint - startpoint)
    co_diff2 = co_diff * (endpoint - midpoint)/(endpoint - startpoint)
    
    t1_idx = np.where(inflow.t < midpoint)[0]
    t1 = inflow.t[t1_idx]
    t2_idx = np.where((inflow.t >= midpoint) & (inflow.t <= endpoint))[0]
    t2_idx = np.insert(t2_idx, 0, t1_idx[-1])
    t2 = inflow.t[t2_idx]
    print(t1_idx, t2_idx)
    
    
    # sinusoid 1
    y1 = co_diff1 * (np.pi / (2 * (midpoint-startpoint)))  * np.sin((np.pi/(midpoint - startpoint)) * t1)
    
    # sinusoid 2
    y2 = co_diff2 * (np.pi / (2 * (endpoint-midpoint)))  * np.sin((np.pi/(endpoint - midpoint)) * (t2 - midpoint))
    print(co_diff1, co_diff2)
    print(np.trapz(y1, t1), np.trapz(y2, t2))
    
    inflow.Q[t1_idx] += y1
    inflow.Q[t2_idx] += y2
    
    
    

    print(len(inflow.Q))
    inflow.update()
    plot_flow(inflow)
    
    print(inflow.mean_inflow * 60/1000)
    
    inflow2 = Inflow(inflow.inflow)
    plot_flow(inflow2)
    
    inflow.write_flow(args.o)
    plot_flow(inflow, save = True, output_file = args.o + '.png')