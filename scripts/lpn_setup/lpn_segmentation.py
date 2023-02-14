# File: lpn_segmentation.py
# File Created: Monday, 13th February 2023 4:03:14 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 14th February 2023 12:09:34 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Creates a LPN

import argparse
import shutil
from pathlib import Path
import subprocess
import os

from svinterface.manager import Manager
from svinterface.core.bc import Inflow, RCR
from svinterface.core.zerod.lpn import LPN

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Creates an LPN")
    parser.add_argument("-i", dest = 'yaml', required = True, help = 'yaml config file')
    parser.add_argument("--f", dest = 'force', action = 'store_true', help = 'force recreation of LPN files')
    args = parser.parse_args()
    
    # manager
    M = Manager(args.yaml)
    
    # construct LPN dir
    lpn_dir = M.root / 'LPN_DIR'
    if args.force and lpn_dir.exists():
        M.unregister("lpn_dir", depth=['workspace'])
        M.unregister("rcrt_file", depth=['workspace'])
        M.unregister("lpn", depth=['workspace'])
        shutil.rmtree(str(lpn_dir))
        # search for rcrt in main dir
        if (lpn_dir / 'rcrt.dat').exists() and not M['options']['tune']:
            M.register("rcrt_file", str(lpn_dir / 'rcrt.dat'), depth=['workspace'])
        
        
    try:
        lpn_dir.mkdir(exist_ok=False)
    except FileExistsError:
        print("New LPN_DIR folder could not be created since one already exists. Use --f flag to overwrite.")
        exit(1)
        
    # register in manager
    M.register(key='lpn_dir',value=str(lpn_dir),depth=['workspace'])
    
    # get files
    files = M['workspace']
    outlet_file = files['outlet_file']
    cent_file = files['centerlines']
    out_dir = str(lpn_dir)
    
    # flowfile, smoothed
    inflow = Inflow.from_file(files['flow_file'], smooth=True, n_points=1000)
    smooth_flow_file = str(lpn_dir / 'inflow_smooth.flow')
    inflow.write_flow(smooth_flow_file)
    
    rcrt_path = lpn_dir / 'rcrt.dat'
    with open(outlet_file, 'r') as ofile:
        outlets = ofile.read().split('\n')[:-1]
    # write rcrt
    if M['options']['tune']:
        # tune mode
        # write an empty if BC does not exist
        bc = RCR()

        for o in outlets:
            bc.add_rcr(o, 0, 0, 0, 0)
        bc.write_rcrt_file(str(lpn_dir))
    else:
        # not tune mode
        if files['rcrt_file'] is None:
            raise FileNotFoundError("rcrt file not found, but BC tuning was not requested.")
        else:
            # copy it into lpn dir
            shutil.copy(files['rcrt_file'], str(rcrt_path))
    # overwrite if one already existed
    M.register('rcrt_file', str(rcrt_path),  ['workspace'])
    
    # get params
    model_name = M['metadata']['model_name']
    # use 1000 timesteps per cycle w/ 6 cycles
    ts = inflow.tc / 1000
    num_ts = 5000

    # run subprocess to generate segmentation
    print("Constructing LPN...", end = '\t', flush = True)
    x = subprocess.run(["simvascular", "--python", "--",  str(Path(__file__).parent / "sv_lpn_segmentation.py"),
                        "-mn", model_name,
                        "-ofc", outlet_file,
                        "-cent", cent_file,
                        "-flow", smooth_flow_file,
                        "-ts", str(ts),
                        "-n_ts", str(num_ts),
                        "-odir", out_dir],
                       capture_output = True)
    
    if x.returncode != 0:
        print("Error")
        # error occured
        print(x)
        exit(1)
    
    print("Done")

    # LPN
    lpn_file = str(lpn_dir / (model_name + '.in'))
    M.register('lpn', lpn_file, ['workspace'])
    
    # Modify LPN
    print("Modifying LPN:")
    lpn = LPN.from_file(lpn_file)
    # add a RCRT map
    print("\tAdding a RCR map...", end = '\t', flush = True)
    lpn.add_rcrt_map(outlets)
    print("Done")
    print("\tConverting to cpp compatible...", end = '\t', flush = True)
    lpn.to_cpp()
    print("Done")
    # update lpn
    lpn.update()
    

    
    # update yaml files at end.
    M.update()
    
    
    