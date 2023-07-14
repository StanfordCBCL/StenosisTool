# File: plot_inflow.py
# File Created: Tuesday, 6th December 2022 5:13:30 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 15th March 2023 3:11:42 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Given Inflow file as argument, plots the inflow.


import sys
from svinterface.plotting import params, plot_flow
from svinterface.core.bc import Inflow



if __name__ == '__main__':
    inflow=Inflow.from_file(sys.argv[1], smooth = False)
    params.set_params()
    plot_flow.plot_flow(inflow, save = False)# save = True, output_file=sys.argv[2])