# File: clear_linear_correction.py
# File Created: Sunday, 4th June 2023 5:26:40 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 9th October 2023 12:11:10 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Clears the linear correction in the LPN_DIR section by overwriting with the previously saved base lpn


from svinterface.core.zerod.lpn import LPN
from svinterface.manager import Manager
import argparse


if __name__ == '__main__':
    
    ## Parser
    parser = argparse.ArgumentParser(description="Clears linear correction")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    args = parser.parse_args()
    
    M = Manager(args.config)
    
    ## Overwrite LPN
    zerod_file = M['workspace']['base_lpn']
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.write_lpn_file(M['workspace']['lpn'])