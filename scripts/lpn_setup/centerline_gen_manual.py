from svinterface.core.polydata import Centerlines
from svinterface.manager import Manager
import argparse
from pathlib import Path

if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description="Generate centerlines")
    parser.add_argument("-mdl", dest = 'mdl', help = 'mdl file')
    parser.add_argument("-vtp", dest = 'vtp', help = 'surface model')
    parser.add_argument("-inlet", dest = 'inlet', help = 'inlet')
    parser.add_argument("-o", dest = 'o', help = 'output file')
    args = parser.parse_args()

    
    # retrieve files
    mdl = args.mdl
    vtp = args.vtp
    inlet = args.inlet
    
    # construct centerlines
    c = Centerlines(None)
    ordered_outlets = c.generate_centerlines(mdl, vtp, inlet, str(args.o))
    
    