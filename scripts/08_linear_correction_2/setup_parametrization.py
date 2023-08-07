

# File: setup_parametrization.py
# File Created: Monday, 17th July 2023 5:41:38 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 17th July 2023 5:53:20 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Sets up a parametrization config section and directory

import argparse
from pathlib import Path

from svinterface.manager import Manager
from svinterface.core.zerod.lpn import LPN

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Set up a parametrization config section and directory")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    parser.add_argument("-s", dest = 'sim', help = 'simulation number to use as base lpn')
    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    # get 0D LPN
    zerod_file = M['simulations'][int(args.sim)]['lpn']
    
    # add param cfg
    if 'parameterization' not in M.yaml or M['parameterization'] is None:
        M.register('parameterization', {})
    
    # add param dir
    param_dir = Path(M['workspace']['root']) / 'parameterization'
    param_dir.mkdir(exist_ok=True)
    M.register('param_dir', str(param_dir), ['workspace'])
    
    # copy base lpn
    base_lpn_path = param_dir / Path(zerod_file).name
    M.register("base_lpn",str(base_lpn_path), ['parameterization'])
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.write_lpn_file(str(base_lpn_path))
    
    # add corrections
    if 'corrections' not in M.yaml['parameterization'] or M.yaml['parameterization']['corrections'] is None:
        M.register('corrections', {}, ['parameterization'])
    
    M.update()