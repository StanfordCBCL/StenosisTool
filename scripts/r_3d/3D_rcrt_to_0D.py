# File: 3d_rcrt_to_1d.py
# File Created: Monday, 31st October 2022 8:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 14th February 2023 1:20:11 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Converts a 3D rcrt file to 0D rcrt file (Only if 3D simulation exists before)
#! Untested.


import re
from svinterface.core.bc import RCR
from pathlib import Path
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(desc = "Converts a 3D rcrt file to 1D rcrt file")
    parser.add_argument("-rcrt", dest = "rcrt_file", help = "3D rcrt file")
    parser.add_argument("-o", dest = "out_dir", help = "output dir to write rcrt.dat")
    parser.add_argument("-svpre", dest = "svpre_file", help = "svpre file for 3D simulation")
    parser.add_argument("-inp", dest = "inp_file", help = "inp file used in 3D simulation")

    args = parser.parse_args()
        
    # read the 3D rcrt
    bc = RCR()
    bc.read_rcrt_file(args.rcrt_file, three_d=True, solver_inp=args.inp_file, svpre=args.svpre_file)
        
    # write the file back out as 0D.
    bc.write_rcrt_file(dirpath = args.out_dir)
    