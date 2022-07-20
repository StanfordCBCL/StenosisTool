import argparse
from functions.centerlines import Centerlines
from functions.utils import *

def dev_centerlines(mod_path, tune_mode = False):
    try:
        files = mod_path.params['files']
        model = mod_path.params['model']
        mdl = files['mdl_file']
        vtp = files['vtp_file']
        inlet = model['inlet']
        
        gen = Centerlines()
        
        if tune_mode:
            if os.path.exists(mod_path.opt_centerlines):
                print('Centerlines file exists already.')
                return False
            
            gen.generate_centerlines(mdl=mdl,
                                vtp=vtp,
                                inlet=inlet,
                                use_entire_tree=False,
                                outlet_names=[model['lpa'], model['rpa']])
        else:
            if os.path.exists(mod_path.model_centerlines):
                print('Centerlines file exists already.')
                return False
        
            gen.generate_centerlines(mdl=mdl,
                                    vtp=vtp,
                                    inlet=inlet,
                                    use_entire_tree=True)
        
        
        if not gen.check_centerlines_data():
            print(mod_path.params['metadata']['name'] + ' failed to generate centerlines.')
            return False
            
        if tune_mode:
            gen.write_centerlines(mod_path.opt_centerlines)
        else:
            gen.write_centerlines(mod_path.model_centerlines)
        return True
    except Exception as e:
        print(e)
        print(mod_path.params['metadata']['name'] + ' failed to generate centerlines.')
        return False

        
        
########
# Main #
########

def tool_main(args):
    raise NotImplementedError()
    # check opt/not opt
    if args.opt and (args.LPA is None or args.RPA is None):
        parser.error('LPA/RPA required with -tune_mode.')
    
        
    
    # tool mode: specified from command line
    mdl = args.mdl_file
    vtp = args.vtp_file
    inlet = args.inlet
    check_exists(mdl)
    check_exists(vtp)
    
    gen = Centerlines()
    
    if args.opt:
        gen.generate_centerlines(mdl=mdl,
                            vtp=vtp,
                            inlet=inlet,
                            use_entire_tree=False,
                            outlet_names=[args.LPA, args.RPA])
    else: 
        gen.generate_centerlines(mdl=mdl,
                            vtp=vtp,
                            inlet=inlet,
                            use_entire_tree=True)
    gen.write_centerlines(args.centerlines_file)

def dev_main(args):
    # dev mode: generates all centerlines
    from functions.project_org import ProjectPath
    
    org = ProjectPath(args.root)
    
    # generate centerlines.
    
    if args.models:
        output = {}
        all_models = org.data_path.models
        for model_name in args.models:
            print('Running for '+ model_name)
            mod_path = org.find_model(model_name)
            if mod_path:
                output[model_name] = dev_centerlines(mod_path, args.tune_mode)
            else:
                output[model_name] = False
    else:
        # run on all models
        output = org.run_all(dev_centerlines, args.tune_mode)
    
    errs = []
    succ = []
    # check which ones have errors
    for mod, val in output.items():
        if val == False:
            errs.append(mod)
        else:
            succ.append(mod)
    print('Models that are successful: ',succ)     
    print('Models with error: ', errs)
    
if __name__ == '__main__':
    
    parser, _, dev, tool = create_parser(description='Constructs centerlines for a model')
    
    
    # dev params
    dev.add_argument('-root', dest = 'root', type = str, default = '.',  help = 'Root to entire project')
    dev.add_argument('-models', dest = 'models', action = 'append', default = [], help = 'Specific models to run')
    dev.add_argument('-tune_mode', dest = 'tune_mode', action = 'store_true', default = False, help = 'flag for whether the centerline generation is for a full model or the optimizer surrogate')

    
    # tool mode
    tool.add_argument('-mdl_file', dest = 'mdl_file', type = str, help = 'File path to mdl file used to map names to face ids')
    tool.add_argument('-vtp_file', dest = 'vtp_file', type = str, help = 'File path to model vtp file')
    tool.add_argument('-inlet', dest = 'inlet', type = str, help = 'inlet face name to not distinguish inlets from outlets')
    tool.add_argument('-tune_mode', dest = 'tune_mode', action = 'store_true', default = False, help = 'flag for whether the centerline generation is for a full model or the optimizer surrogate')
    tool.add_argument('-LPA', dest = 'lpa', default = None, help = 'LPA required for tune_mode')
    tool.add_argument('-RPA', dest = 'rpa', default = None, help = 'RPA required for tune_mode')
    tool.add_argument('-output_file', dest = 'centerlines_file', type = str, help = 'Destination vtp file to write centerlines to. Include vtp')
    
    args = parser.parse_args()
    
    
        
    if args.mode == 'tool':
        tool_main(args)
        
    elif args.mode == 'dev':
        dev_main(args)
        
        
    
        
        
        
        
        
    
    
    