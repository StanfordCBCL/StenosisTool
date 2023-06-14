# File: scale_flow_amp.py
# File Created: Wednesday, 31st May 2023 2:20:24 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 31st May 2023 2:56:58 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Scales flow to a certain systolic flow without changing diastolic flow. Does so by multiplying all values by a ratio r s.t. after shifting the inflow down by the Q_dia * (r - 1), we still achieve the target systolic flow


from svinterface.core.bc import Inflow
from svinterface.plotting.plot_flow import plot_flow

import argparse
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Add sine to plot")
    
    parser.add_argument("-i", help = 'inflow waveform')
    parser.add_argument("-o", help = "inflow waveform scaled out file")
    parser.add_argument("--s", dest = "sys", type=float, help = "systolic flow target (mL/s)")
    
    args = parser.parse_args()

    inflow = Inflow.from_file(args.i, smooth=False)
    # plot_flow(inflow)
    
    print("Old min:", inflow.min_inflow)
    print("Old max:", inflow.max_inflow)
    print("Old mean:", inflow.mean_inflow)
    # compute ratio r
    x = (args.sys - inflow.max_inflow) * inflow.min_inflow / (inflow.max_inflow - inflow.min_inflow)
    r = (args.sys + x) / inflow.max_inflow
    
    # multiply by r
    inflow.Q *= r
    # shift down by Q_dia * (r - 1)
    inflow.Q -= (inflow.min_inflow * (r - 1))
    
    
    inflow.update()
    print("New min:", inflow.min_inflow)
    print("New max:", inflow.max_inflow)
    print("New mean:", inflow.mean_inflow)
    
    plot_flow(inflow, save = True, output_file = args.o + '.png')
    inflow.write_flow(args.o)
    