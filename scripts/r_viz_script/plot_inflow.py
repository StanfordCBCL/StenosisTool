# File: plot_inflow.py
# File Created: Tuesday, 6th December 2022 5:13:30 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 26th January 2023 4:19:56 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Given Inflow file as argument, plots the inflow.


import sys
from sgt.core.flow import Inflow


if __name__ == '__main__':
    inflow=Inflow.from_file(sys.argv[1], smooth = False)
    inflow.plot_flow(save = False)