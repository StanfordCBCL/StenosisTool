# File: sv_dev_centerline_gen.py
# File Created: Thursday, 28th July 2022 3:28:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 17th October 2022 2:54:59 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Using Simvascular VMTK, generate centerlines for a particular 3D geometry


import sys
import os
# append the src path, since it uses Simvascular's python.
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.centerlines import Centerlines
from src.misc import create_dev_parser
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
            if os.path.exists(model.model_centerlines):
                print('Model centerlines file exists already.')
                return False
    
        gen.generate_centerlines(mdl=mdl,
                                vtp=vtp,
                                inlet=inlet,
                                use_entire_tree=True)
        
        
        if not gen.check_centerlines_data():
            print(model.info['metadata']['name'] + ' failed to generate centerlines.')
            return False
            
        gen.write_centerlines(model.model_centerlines)
        return True
    except Exception as e:
        print(e)
        print(model.info['metadata']['name'] + ' failed to generate centerlines.')
        return False


def determine_success(output):
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
    else:
        model_list = list(org.model_names)
    
    # iterate through desired models
    for model_name in model_list:
        model = org.find_model(model_name)
        # if model does not exist, skip
        if model is None:
            output[model_name] = False
            continue
    
        output[model_name] = construct_full_centerlines(model, force=args.force)

    err, succ = determine_success(output)
    
    print('Models that are successful: ',succ)     
    print('Models with error: ', err)

    
if __name__ == '__main__':
    
    dev_parser = create_dev_parser(desc='Constructs centerlines for a model')
    
    # dev params
    dev_parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'whether to force overwriting existing centerlines')

    args = dev_parser.parse_args()
    main(args)
        
        
    
        
        
        
        
        
    
    
    