
from src.misc import create_dev_parser
from src.solver import Solver0D
from src.bc import BoundaryConditions
from src.data_org import DataPath
from src.file_io import parse_mdl, write_json, check_many

import os
import shutil
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

########
# Main #
########

def main(args):
    
    org = DataPath(args.root)
    
    if args.models:
        model_list = args.models
    else:
        model_list = list(org.model_names)
    
    for model_name in model_list:
        print(f'Configuring base solver for {model_name}...', end = '\t', flush=True)
        model = org.find_model(model_name)
        
        # func to join to base solver
        join = lambda file: os.path.join(model.base_solver_dir, file)
        
        # Convert Model to cpp 
        #! May be able to remove with updates to C solver
        solver = Solver0D()
        solver.read_solver_file(model.base_model_solver)
        solver.to_cpp()
        
        # add a BC map
        bc = BoundaryConditions()
        bc.read_rcrt_file(join('rcrt.dat'))
        solver.write_bc_map(bc)
        
        # write model units and name
        solver.simulation_params['units'] = model.info['model']['units']
        solver.simulation_params['age'] = int(model.info['metadata']['age'])
        solver.simulation_params['gender'] = model.info['metadata']['gender']
        
        # plot inflow
        solver.inflow.plot_flow(join('inflow.png'))
        
        # write cpp solver
        solver.write_solver_file(model.base_model_solver)
        
        # write inlet/outlet mappings file
        outlet_mappings = parse_mdl(model.info['files']['mdl_file'])
        inlet_mapping = {model.info['model']['inlet']: outlet_mappings[model.info['model']['inlet']]}
        del outlet_mappings[model.info['model']['inlet']]
        write_json(join('inlet_mapping.dat'), inlet_mapping)
        write_json(join('outlet_mapping.dat'), outlet_mappings)
        
        # copy centerlines
        shutil.copy(model.model_centerlines, join('model_centerlines.vtp'))
        
        # if stenosis model, write the rcrt into 0D format
        if model.type == 'stenosis':
            try:
                convert_old_rcrt(old_rcrt_file=model.info['files']['rcrt3d_file'],
                                solver3d_file=model.info['files']['solver3d_file'],
                                pre_file=model.info['files']['pre3d_file'],
                                out_dir=model.base_solver_dir,
                                prefix=model.info['model']['prefix'])
            except Exception as e:
                print(f'Failed: {e}')
                continue
        elif model.type == 'healthy':
            
            # copy CapInfo
            if not os.path.exists(join('CapInfo')):
                shutil.copy(model.info['files']['cap_info'], join('CapInfo'))
            
        print('Done')
        
        
        
        
        

if __name__ == '__main__':
    
    dev_parser = create_dev_parser(desc = 'set up base solver with all necessary files for other non-dev scripts')
    
    args = dev_parser.parse_args()
    main(args)
    