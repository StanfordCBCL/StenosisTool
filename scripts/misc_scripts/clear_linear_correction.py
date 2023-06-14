from svinterface.core.zerod.lpn import LPN
from svinterface.manager import Manager
from svinterface.core.bc import RCR
import argparse


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Clears linear correction")
    parser.add_argument("-i", dest = 'config', help = 'config.yaml file')
    
    args = parser.parse_args()
    
    
    M = Manager(args.config)
    
    zerod_file = M['workspace']['lpn']
    
    # get LPN and convert to BVJ
    zerod_lpn = LPN.from_file(zerod_file)
    zerod_lpn.to_cpp(normal = False) # resets to all 0's
    # reloads the rcrts from previously tuned
    rcr = RCR()
    rcr.read_rcrt_file(M['workspace']['rcrt_file'])
    zerod_lpn.update_rcrs(rcr)
    
    zerod_lpn.update()