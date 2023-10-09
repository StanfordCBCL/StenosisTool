# File: invert_flow.py
# File Created: Wednesday, 31st May 2023 5:14:36 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 9th October 2023 12:12:42 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Inverts a given flow.

from svinterface.core.bc import Inflow
import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Invert flow")
    parser.add_argument("-i", help = 'inflow waveform')
    parser.add_argument("-o", help = "outflow file of inverted inflow")
    args = parser.parse_args()
    
    ## Load in inverse flow
    i = Inflow.from_file(args.i, inverse = True, smooth = False)
    i.write_flow(args.o)
    print(f"Inverted inflow written to {args.o}")