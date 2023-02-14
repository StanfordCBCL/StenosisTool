# File: setup_workspace.py
# File Created: Monday, 31st October 2022 6:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 13th February 2023 9:38:55 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Sets up a results workspace from a SV Model Directory

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
    
    parser = argparse.ArgumentParser(description = 'Sets up workspace from SV Model Directory')
    
    parser.add_argument("-i", dest = 'yaml', required = True, help = "SV Yaml file")
    parser.add_argument("-o", dest = 'outdir', required = True, help = "Output directory to construct workspace in.")
    parser.add_argument('--f', dest = 'force', action = 'store_true', default = False, help = 'force re-set up workspace')
    args = parser.parse_args()
    
    # Manager
    sv = Manager(args.yaml)
    
    # outdir path
    outdir = Path(args.outdir)
    
    # if force is given, del old outdir if it exists
    if args.force and outdir.exists():
        shutil.rmtree(str(outdir))
        
    # construct dir, should not already exist
    try:
        outdir.mkdir(parents=True, exist_ok = False)
    except FileExistsError as e:
        print("Workspace already exists. To overwrite, use --f flag.")
        exit(1)
        
    # add from workspace to file
    print("Copying workspace files...", end = '\t', flush = True)

    # update yaml
    assert sv['metadata']['model_name'] != None, "Please add a model name to your yaml file."
    model_name = sv['metadata']['model_name']

    # copy stuff over
    workspace = sv['workspace']
    
    # surface model
    copy_item(sv, workspace, key = 'surface_model', filename = model_name + '.vtp')
    # mdl
    copy_item(sv, workspace, key = 'mdl', filename = model_name + '.mdl') 
    # flow
    copy_item(sv, workspace, key = 'flow_file', filename = 'inflow.flow')
    # BC if precomputed
    if workspace['rcrt_file'] != None:
        copy_item(sv, workspace, key = 'rcrt_file', filename = 'rcrt.dat')
    # if Capinfo is found
    if workspace['capinfo']:
        copy_item(sv, workspace, key = 'capinfo', filename = 'CapInfo')
        
        
    # write root
    workspace['root'] = str(outdir)
    
    print("Done")
    sv.write(str(outdir/'config.yaml'))
    print('Successfully constructed workspace:', args.outdir)
    
    