# File: centerline_gen_manual.py
# File Created: Tuesday, 13th June 2023 6:07:27 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 13th July 2023 1:50:48 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Used to generate centerlines for any SV Model. Does not provide the order of the outlets.

from svinterface.core.polydata import Centerlines
import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Generate centerlines")
    parser.add_argument("-mdl", dest='mdl', required=True, help='Path to mdl file associated w/ model')
    parser.add_argument("-vtp", dest='vtp', required=True, help='Path to vtp surface model file')
    parser.add_argument("-inlet", dest='inlet', required=True, help='Name of the inlet cap')
    parser.add_argument("-o", dest = 'outdir', help = 'Path to output the file.')
    args = parser.parse_args()
    
    # construct centerlines
    c = Centerlines(None)
    ordered_outlets = c.generate_centerlines(args.mdl, args.vtp, args.inlet, args.output)
    
    