# File: setup_workspace.py
# File Created: Monday, 31st October 2022 6:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 14th September 2023 2:51:50 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Sets up a results workspace from a SV Model Directory. Requires a config file from generate_sv_config.py

from svinterface.manager import Manager

from pathlib import Path
import shutil
import argparse

def copy_item(sv, workspace, key, filename):
    """Copies a particular item over
    """
    sm = outdir / filename
    shutil.copy(str(sv.root / workspace[key]), str(sm))
    workspace[key] = str(sm)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Sets up a results workspace from a SV Model Directory. Requires a config file from generate_sv_config.py')
    parser.add_argument("-i", dest='config', required=True, help="Config file from generate_sv_config.py")
    parser.add_argument("-o", dest='outdir', required=True, help="Path to the directory to construct workspace in. If it already exists, an error is thrown.")
    args = parser.parse_args()
    
    sv = Manager(args.config)
    
    ## Construct Directory
    outdir = Path(args.outdir)
    try:
        outdir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as e:
        print("Workspace already exists, please remove it first.")
        exit(1)
        
    ## Copy over files
    print("Copying workspace files...", end='\t', flush=True)

    assert sv['metadata']['model_name'] != None, "Please add a model name to your yaml file."
    model_name = sv['metadata']['model_name']
    workspace = sv['workspace']

    copy_item(sv, workspace, key='surface_model', filename=model_name+'.vtp')
    copy_item(sv, workspace, key='mdl', filename=model_name+'.mdl') 
    copy_item(sv, workspace, key='flow_file', filename='inflow.flow')
    if workspace['rcrt_file'] != None:
        copy_item(sv, workspace, key='rcrt_file', filename='rcrt.dat')
    if workspace['capinfo']:
        copy_item(sv, workspace, key='capinfo', filename='CapInfo')
    
    workspace['root'] = str(outdir)
    print("Done")
    
    ## Write updated config
    sv.write(str(outdir/'config.yaml'))
    print('Successfully constructed workspace:', args.outdir)
    
    