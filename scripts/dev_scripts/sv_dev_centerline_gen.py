# File: sv_dev_centerline_gen.py
# File Created: Thursday, 28th July 2022 3:28:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 19th October 2022 7:31:27 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Using Simvascular VMTK, generate centerlines for a particular 3D geometry


import sys
import os
from pathlib import Path
# append the src path, since it uses Simvascular's python.
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.polydata import Centerlines
from src.parser import DevParser
from src.data_org import DataPath, ModelPath

def construct_full_centerlines(model: ModelPath, force = False):
    ''' Constructs full centerlines of model'''
    try:
        files = model.info['files']
        m = model.info['model']
        mdl = files['mdl_file']
        vtp = files['vtp_file']
        inlet = m['inlet']
        
        gen = Centerlines()
        
        if not force:
            if model.model_centerlines.exists():
                print('Model centerlines file exists already.')
                return False
    
        gen.generate_centerlines(mdl=str(mdl),
                                vtp=str(vtp),
                                inlet=str(inlet),
                                use_entire_tree=True)
        
        
        if not gen.check_centerlines_data():
            print(model.model_name + ' failed to generate centerlines.')
            return False
            
        gen.write_polydata(str(model.model_centerlines))
        return True
    
    except Exception as e:
        print(e)
        print(model.model_name + ' failed to generate centerlines.')
        return False


def determine_success(output):
    ''' check which models succeeded or failed
    '''
    errs = []
    succ = []
    # check which ones have errors
    for mod, val in output.items():
        if val == False:
            errs.append(mod)
        else:
            succ.append(mod)
    return errs, succ

########
# Main #
########

def main(args):
    
    org = DataPath(args.root)

    output = {}

    # get model list
    if args.models:
        model_list = args.models
    # otherwise all models
    else:
        model_list = sorted(org.model_names)
    
    # iterate through desired models
    for model_name in model_list:
        model = org.find_model(model_name)
        # if model does not exist, skip
        if model is None:
            output[model_name] = False
            continue
        
        # construct centerlines
        output[model_name] = construct_full_centerlines(model, force=args.force)

    # report success/failure
    err, succ = determine_success(output)
    print('Models that are successful: ',succ)     
    print('Models with error: ', err)

    
if __name__ == '__main__':
    
    dev_parser = DevParser(desc='Constructs centerlines for a model')
    dev_parser.parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'whether to force overwriting existing centerlines')

    args = dev_parser.parse_args()
    main(args)
        
        
    
        
        
        
        
        
    
    
    