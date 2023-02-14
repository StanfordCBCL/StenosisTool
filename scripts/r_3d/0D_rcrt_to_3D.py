# File: 0D_rcrt_to_3D.py
# File Created: Monday, 31st October 2022 8:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 14th February 2023 12:29:22 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Converts a 0D rcrt file to 3D rcrt file.


import re
from svinterface.core.bc import RCR
from pathlib import Path
import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = "Converts a 3D rcrt file to 1D rcrt file")
    parser.add_argument("-rcrt", dest = "rcrt_file", help = "0D rcrt file")
    parser.add_argument("-o", dest = "out_dir", help = "output dir to write rcrt.dat")
    parser.add_argument("-svpre", dest = "svpre_file", help = "svpre file for 3D simulation")
    parser.add_argument("-inp", dest = "inp_file", help = "inp file used in 3D simulation")

    args = parser.parse_args()
    
    bc = RCR()
    
    bc.read_rcrt_file(args.rcrt_file)
    bc.sort_for_3d(args.inp_file, args.svpre_file)
    bc.write_rcrt_file(args.out_dir, as_3d=True)