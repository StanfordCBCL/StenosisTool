# File: sv_dev_centerline_gen.py
# File Created: Thursday, 28th July 2022 3:28:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 11th November 2022 12:38:35 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Using Simvascular VMTK, generate centerlines for a particular 3D geometry


import sys
import os
from pathlib import Path
# append the src path, since it uses Simvascular's python.
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from sgt.core.polydata import Centerlines
from sgt.utils.parser import Parser


def construct_full_centerlines(mdl, vtp, inlet, out_file):
    
    ''' Constructs full centerlines of model'''
    try:
        gen = Centerlines()
        gen.generate_centerlines(mdl=str(mdl),
                                vtp=str(vtp),
                                inlet=str(inlet),
                                use_entire_tree=True)
        
        
        if not gen.check_centerlines_data():
            print(Path(vtp).stem + ' failed to generate centerlines: Centerlines data incomplete.')
            return False
            
        gen.write_polydata(out_file)
        return True
    
    except Exception as e:
        print(e)
        print(Path(vtp).stem + ' failed to generate centerlines.')
        return False


    
if __name__ == '__main__':
    
    parser = Parser(desc='Constructs centerlines for a model')
    
    parser.parser.add_argument('-mdl', help = 'mdl file path')
    parser.parser.add_argument('-vtp', help = 'surface model to construct centerlines from')
    parser.parser.add_argument('-inlet', help = 'inlet cap')
    parser.parser.add_argument('-o', help = 'output file')

    args = parser.parse_args()
    if construct_full_centerlines(args.mdl, args.vtp, args.inlet, args.o):
        print("Centerline generation successful")
        
    else:
        print("Centerline generation failed.")
    
        
        
    
        
        
        
        
        
    
    
    