# File: 0D_rcrt_to_3D.py
# File Created: Monday, 31st October 2022 8:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 7:18:45 pm
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
    parser.add_argument("-o", dest = "out_dir", help = "output dir to write rcrt.dat (must contain an svpre and an inp file)")
    parser.add_argument('-rep',required=False, default = [], action='append',nargs=2 ,help='For stented models, certain names may be replaced with a different one.  old name, new name')
    args = parser.parse_args()
    
    ## Load Boundary Conditions
    bc = RCR()
    bc.read_rcrt_file(args.rcrt_file)
    
    # find inp and svpre files
    inp_file = None
    svpre_file = None
    out_dir = Path(args.out_dir)
    for files in out_dir.iterdir():
        if files.suffix == '.inp':
            inp_file = str(files)
        elif files.suffix == '.svpre':
            svpre_file = str(files)
    
    if inp_file is None:
        print(".inp file not found. Error")
        exit(1)
    elif svpre_file is None:
        print(".svpre file not found. Error")
        exit(1)
    
    # if replacing name is necessary for diseased simulation
    for old, new in args.rep:
        bc.bc_list[new] = bc.bc_list[old]
        bc.bc_list[new]['faceID'] = new
        del bc.bc_list[old]
    
    # Sort for 3D and write file
    bc.sort_for_3d(inp_file, svpre_file)
    bc.write_rcrt_file(args.out_dir, as_3d=True)