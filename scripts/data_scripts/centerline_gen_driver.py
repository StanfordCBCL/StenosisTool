
#! SINCE THIS NEEDS TO BE RUN WITH SIMVASCULAR PYTHON, THIS HACK IS NECESSARY UNTIL A MORE ELEGANT SOLUTION IS FOUND
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.centerlines import Centerlines
from src.misc import create_parser
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

def tool_main(args):
    raise NotImplementedError
    # TODO: Take in model root path, construct a ModelPath Object
    # TODO: Feed it into construct_tuning or construct full depending on which flag
    

def dev_main(args):
    
    org = DataPath(args.root)
    # generate centerlines.
    print(org)
    
    output = {}
    for model_name in args.models:
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
    
    parser, dev, tool = create_parser(desc='Constructs centerlines for a model')
    
    # dev params
    dev.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    dev.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'whether to force overwriting existing centerlines')

    args = parser.parse_args()
        
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)
        
        
    
        
        
        
        
        
    
    
    