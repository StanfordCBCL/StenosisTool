# File: scale_flow_sine.py
# File Created: Monday, 22nd May 2023 12:02:17 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 9th October 2023 12:15:00 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Scales inflow to be appropriate using 2 sinusoids. Does not shift systolic

from svinterface.core.bc import Inflow
from svinterface.plotting.plot_flow import plot_flow
import argparse
import numpy as np

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Add sine to plot")
    parser.add_argument("-i", help = 'inflow waveform')
    parser.add_argument("-o", help = "inflow waveform scaled out file")
    parser.add_argument("--co", dest = "co", type=float, help = "cardiac outflow")
    args = parser.parse_args()

    inflow = Inflow.from_file(args.i, smooth=False)
    
    ## Get points for Sinusoids
    startpoint = 0
    midpoint = inflow.max_inflow_t
    endpoint = inflow.min_inflow_t
    
    ## Compute CO Diff
    # print(args.co * 1000/60, inflow.mean_inflow)
    co_diff = args.co * 1000/60 - inflow.mean_inflow
    
    ## Split Ratio by ranges
    co_diff1 = co_diff * (midpoint - startpoint)/(endpoint - startpoint)
    co_diff2 = co_diff * (endpoint - midpoint)/(endpoint - startpoint)
    
    ## Compute indices
    t1_idx = np.where(inflow.t < midpoint)[0]
    t1 = inflow.t[t1_idx]
    t2_idx = np.where((inflow.t >= midpoint) & (inflow.t <= endpoint))[0]
    t2_idx = np.insert(t2_idx, 0, t1_idx[-1])
    t2 = inflow.t[t2_idx]
    # print(t1_idx, t2_idx)
    
    
    # Sinusoid 1
    y1 = co_diff1 * (np.pi / (2 * (midpoint-startpoint)))  * np.sin((np.pi/(midpoint - startpoint)) * t1)
    
    # Sinusoid 2
    y2 = co_diff2 * (np.pi / (2 * (endpoint-midpoint)))  * np.sin((np.pi/(endpoint - midpoint)) * (t2 - midpoint))
    # print(co_diff1, co_diff2)
    # print(np.trapz(y1, t1), np.trapz(y2, t2))

    # Add Sinusoids
    inflow.Q[t1_idx] += y1
    inflow.Q[t2_idx] += y2
    
    ## Plotting and Saving
    inflow.update()
    plot_flow(inflow)
    # print(inflow.mean_inflow * 60/1000)
    inflow2 = Inflow(inflow.inflow)
    plot_flow(inflow2)
    inflow.write_flow(args.o)
    plot_flow(inflow, save = True, output_file = args.o + '.png')