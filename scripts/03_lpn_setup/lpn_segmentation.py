# File: lpn_segmentation.py
# File Created: Monday, 13th February 2023 4:03:14 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 13th February 2024 2:25:08 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Constructs an LPN from the recorded centerlines in the manager. Uses preset parameters only.

import argparse
import shutil
from pathlib import Path
import subprocess

from svinterface.manager import Manager
from svinterface.core.bc import Inflow, RCR
from svinterface.core.zerod.lpn import LPN

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Constructs an LPN from the recorded centerlines in the manager. Uses preset parameters only.")
    parser.add_argument("-i", dest = 'config', required = True, help = 'Config file')
    parser.add_argument("--f", dest = 'force', action = 'store_true', help = 'force recreation of LPN files')
    args = parser.parse_args()
    
    ## Manager
    M = Manager(args.config)
    
    ## Construct LPN dir
    lpn_dir = M.root / 'LPN_DIR'
    # reset if forced to.
    if args.force and lpn_dir.exists():
        M.unregister("lpn_dir", depth=['workspace'])
        M.unregister("rcrt_file", depth=['workspace'])
        M.unregister("lpn", depth=['workspace'])
        shutil.rmtree(str(lpn_dir))
        # search for rcrt in main dir
        if (lpn_dir / 'rcrt.dat').exists() and not M['options']['tune']:
            M.register("rcrt_file", str(lpn_dir / 'rcrt.dat'), depth=['workspace'])
    # create
    try:
        lpn_dir.mkdir(exist_ok=False)
    except FileExistsError:
        print("New LPN_DIR folder could not be created since one already exists. Use --f flag to overwrite.")
        exit(1)
        
    # register in manager
    M.register(key='lpn_dir',value=str(lpn_dir),depth=['workspace'])
    
    
    ## LPN Construction
    # get files
    files = M['workspace']
    outlet_file = files['outlet_file']
    cent_file = files['centerlines']
    out_dir = str(lpn_dir)
    
    # Retrieve inflow file, but smooth it.
    inflow = Inflow.from_file(files['flow_file'], smooth=True, n_points=1000)
    smooth_flow_file = str(lpn_dir / 'inflow_smooth.flow')
    inflow.write_flow(smooth_flow_file)
    
    # Create RCRT file
    rcrt_path = lpn_dir / 'rcrt.dat'
    with open(outlet_file, 'r') as ofile:
        outlets = ofile.read().split('\n')[:-1]
    # Write rcrt
    if M['options']['tune']:
        # if tuning is required, write an empty RCRT file
        bc = RCR()
        for o in outlets:
            bc.add_rcr(o, 0, 0, 0, 0)
        bc.write_rcrt_file(str(lpn_dir))
    else:
        # If tuning is not required, search for existing RCRT file.
        if files['rcrt_file'] is None:
            raise FileNotFoundError("rcrt file not found, but BC tuning was not requested.")
        else:
            shutil.copy(files['rcrt_file'], str(rcrt_path))
    M.register('rcrt_file', str(rcrt_path),  ['workspace'])
    
    # get parameters
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

    ## LPN register.
    lpn_file = str(lpn_dir / (model_name + '.in'))
    M.register('lpn', lpn_file, ['workspace'])
    
    ## Modify constructed LPN for future work
    print("Modifying LPN:")
    lpn = LPN.from_file(lpn_file)
    
    # add a RCRT map
    print("\tAdding a RCR map...", end = '\t', flush = True)
    lpn.add_rcrt_map(outlets)
    print("Done")
    # convert to a CPP compatible version
    print("\tConverting to cpp compatible...", end = '\t', flush = True)
    lpn.to_cpp()
    print("Done")

    lpn.update()
    
    ## Append the rcrs and register the base lpn if tuning is not required
    if not M['options']['tune']:
        bcs = RCR()
        bcs.read_rcrt_file(str(rcrt_path))
        lpn.update_rcrs(bcs)
        
        base_lpn_path = Path(M['workspace']['root']) / 'base_lpn.in'
        M.register('base_lpn', str(base_lpn_path), depth = ['workspace'])
        lpn.write_lpn_file(str(base_lpn_path))
    
    # update Manager
    M.update()
    
    
    