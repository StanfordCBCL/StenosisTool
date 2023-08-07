from svinterface.core.zerod.lpn import LPN
from svinterface.manager import Manager
from svinterface.core.bc import RCR
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Clears linear correction")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    
    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    # overwrite base to lpn
    zerod_file = M['workspace']['base_lpn']
    
    # get lpn
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.write_lpn_file(M['workspace']['lpn'])