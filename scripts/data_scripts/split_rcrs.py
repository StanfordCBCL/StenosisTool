from src.misc import create_parser
from src.data_org import DataPath
from src.bc import BoundaryConditions

import numpy as np
import os


#############
# Split RCR #
#############

def load_area_file(area_filepath):
    with open(area_filepath, 'r') as afile:
        areas = {}
        afile.readline() # ignore first comment line
        for line in afile:
            line = line.rstrip().split()
            areas[line[0]] = float(line[1])
    return areas

def total_area(areas):
    return sum(list(areas.values()))

def validate_caps(areas):
    ''' confirm that all caps have either lpa or rpa in them'''
    names = list(areas.keys())
    for name in names:
        if 'lpa' not in name.lower() and 'rpa' not in name.lower():
            raise ValueError('Unable to identify RPA vs. LPA caps: please rename caps to contain lpa or rpa')
    return 

def split_rpa_lpa(areas):
    ''' splits areas between lpa and rpa areas '''
    lpa_areas = {}
    rpa_areas = {}
    
    validate_caps(areas)
    
    for name, area in areas.items():
        if 'lpa' in name.lower():
            lpa_areas[name] = area
        elif 'rpa' in name.lower():
            rpa_areas[name] = area
    return lpa_areas, rpa_areas


def split_bc(areas,  x):
    
    def Rpi(Ai, A, Rp):
        return (A / Ai) * Rp
    
    def Rdi(Ai, A, Rd):
        return (A / Ai) * Rd
    
    def Ci(Ai, A, C):
        return (Ai/A) * C
    
    lpa_areas, rpa_areas = split_rpa_lpa(areas)
    
    lpa_A = total_area(lpa_areas)
    rpa_A = total_area(rpa_areas)
    
    rcrs = {}
    
    
    for name, area in lpa_areas.items():
        rcrs[name] = {'Rp': Rpi(area, lpa_A, x[0]),
                     'C': Ci(area, lpa_A, x[1]),
                     'Rd': Rdi(area, lpa_A, x[2]) }
    
    for name, area in rpa_areas.items():
        rcrs[name] = {'Rp': Rpi(area, rpa_A, x[3]),
                     'C': Ci(area, rpa_A, x[4]),
                     'Rd': Rdi(area, rpa_A, x[5]) }
    
    return rcrs

def calc_rcrs(results, inlet, area_file, out_dir):
    ''' Splits RCRS '''
    
    areas = load_area_file(area_file)
    
    del areas[inlet]
    
    rcrs = split_bc(areas, results['x'])
    
    bcs = BoundaryConditions()
    
    for name, vals in rcrs.items():
        bcs.add_rcr(face_name=name, Rp = vals['Rp'], C = vals['C'], Rd = vals['Rd'], Pd = results['Pd'])
    
    bcs.write_rcrt_file(out_dir)
    
    return rcrs

#########
# Mains #
#########

def tool_main(args):
    # TODO: CREATE MODEL AND SPLIT
    raise NotImplementedError

def dev_main(args):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        model = org.find_model(model_name)
        
        if 'cap_info' in model.info['files']:
            area_file = model.info['files']['cap_info']
        else:
            print('Missing CapInfo for ' + model_name)
            continue
        results_file = os.path.join(os.path.dirname(model.tune_solver), model.info['metadata']['name'] + '_optimize_results.npy')
        results = np.load(results_file, allow_pickle=True ).item()
        
        calc_rcrs(results=results, 
                  inlet=model.info['model']['inlet'],
                  area_file=area_file,
                  out_dir=model.solver_dir)
            


if __name__ == '__main__':
    
    parser, dev, tool = create_parser(desc = 'Splits optimized results in rcrt.dat')
    
    dev.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    
    args = parser.parse_args()
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)