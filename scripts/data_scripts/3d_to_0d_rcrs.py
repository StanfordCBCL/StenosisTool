
from lib2to3.pytree import convert
from src.misc import create_parser
from src.data_org import DataPath
from src.bc import BoundaryConditions
import os
import re

#############
# Functions #
#############

def read_pre_file(pre_file):
    
    mapping = {}
    with open(pre_file, 'r') as pfile:
        for line in pfile:
            if line.startswith('set_surface_id_vtp'):
                line = line.split()
                face_id = int(line[2])
                face_name = os.path.splitext(os.path.basename(line[1]))[0]
                mapping[face_id] = face_name
    return mapping

def read_solver_file(solver3d_file):
    with open(solver3d_file, 'r') as sfile:
        for line in sfile:
            if re.search('List of RCR Surfaces:', line):
                used_values = line.split(':')[1].split()
                break

    return [int(i) for i in used_values]

    
                
def convert_old_rcrt(old_rcrt_file, solver3d_file, pre_file, out_dir, prefix = '_'):
    
    face_mapping = read_pre_file(pre_file)
    used_surfaces = read_solver_file(solver3d_file)

    
    with open(old_rcrt_file, 'r') as rfile:
        keyword = rfile.readline()
        cur_cap_idx = 0
        bcs = BoundaryConditions()
        while True:
            #print(cur_cap_idx, face_mappings[ids[cur_cap_idx]])
            
            
            tmp = rfile.readline()
            if tmp == keyword:
                face_name = prefix + face_mapping[used_surfaces[cur_cap_idx]]
                Rp = float(rfile.readline())
                C = float(rfile.readline())
                Rd = float(rfile.readline())
                p0 = float(rfile.readline().strip().split()[1])
                p1 = float(rfile.readline().strip().split()[1])
                assert p0 == p1, 'Cannot handle time-dependent reference pressure'
                Pd = (float(p1))
                
                # add rcr bc
                bcs.add_rcr(face_name = face_name, Rp = Rp, C = C, Rd = Rd, Pd = Pd)
            
            if cur_cap_idx > len(used_surfaces):
                raise RuntimeError('More rcrt BCs than possible caps')
            
            cur_cap_idx += 1
            if len(tmp) == 0:
                break
            
        if cur_cap_idx <= len(used_surfaces):
            raise RuntimeError('Insuffienct rcrt BCs to match caps')
                
    bcs.write_rcrt_file(out_dir)



#########
# Mains #
#########
def tool_main(args):
    raise NotImplementedError

def dev_main(args):
    
    org = DataPath(args.root)
    
    for model_name in args.models:
        model = org.find_model(model_name)
        
        fail = False
        # check needed params are there
        for filepath in ['solver3d_file', 'pre3d_file', 'rcrt3d_file']:
            if filepath not in model.info['files']:
                print(f'{filepath} not found')
                fail = True
        if fail:
            print('Failure to find files: skipping')
            continue
        
        convert_old_rcrt(model.info['files']['rcrt3d_file'],
                         model.info['files']['solver3d_file'],
                         model.info['files']['pre3d_file'],
                         model.solver_dir,
                         args.prefix + '_')
        
        

if __name__ == '__main__':
    
    parser, dev, tool = create_parser(desc = 'converts 3d rcrt to 0d')
    
    dev.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    dev.add_argument('-prefix', default = '', help = 'prefix to add before bc name')
    
    args = parser.parse_args()
    
    if args.mode == 'tool':
        tool_main(args)
    elif args.mode == 'dev':
        dev_main(args)