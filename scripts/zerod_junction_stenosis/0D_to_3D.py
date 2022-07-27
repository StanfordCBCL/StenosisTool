from src.bc import BoundaryConditions
import sys
import os
import re
import shutil

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

def write_solver_file(solver3d_file, rcr_surfaces):
    with open(solver3d_file, 'r') as sfile:
        file = ''
        for line in sfile:
            
            if re.search('Number of Resistance Surfaces:', line):
                file += 'Number of RCR Surfaces: ' + str(len(rcr_surfaces)) + '\n'
            elif re.search('List of Resistance Surfaces', line):
                file += 'List of RCR Surfaces: ' + ' '.join([str(x) for x in rcr_surfaces]) + '\n'
            elif re.search('Resistance Values', line):
                file += 'RCR Values From File: True\n'
            elif re.search('Number of Timesteps between Restarts:', line):
                file += 'Number of Timesteps between Restarts: 200\n'
            else:
                file += line
        
    with open(solver3d_file,'w') as sfile:
        sfile.write(file)

if __name__ == '__main__':
    
    rcrt_file = sys.argv[1]
    threed_dir = sys.argv[2]
    
    solver_file = os.path.join(threed_dir, 'solver.inp')
    backup = os.path.splitext(solver_file)[0] + '_backup.inp'
    if not os.path.exists(backup):
        shutil.copy(solver_file, backup)
    for file in os.listdir(threed_dir):
        if os.path.splitext(file)[-1] == '.svpre':
            
            svpre_file = os.path.join(threed_dir, file)

    bc = BoundaryConditions()
    bc.read_rcrt_file(rcrt_file)
    
    bc_map = {}
    
    for bc_val in bc.bc_list:
        bc_map[bc_val['faceID']] = bc_val
        
    unused_ids = []
    face_mapping = read_pre_file(svpre_file)
    for faceid, bcname in face_mapping.items():
        if bcname not in bc_map:
            unused_ids.append(faceid)
    for faceid in unused_ids:
        del face_mapping[faceid]
    
    rcr_surfaces = sorted(list(face_mapping.keys()))
    # write solver file
    write_solver_file(solver_file, rcr_surfaces)
    new_rcrt = os.path.join(threed_dir, 'rcrt.dat')
    # write rcrt
    with open(new_rcrt, 'w') as rcrt_file:
        rcrt_file.write('2\n')
        for rcr_id in rcr_surfaces:
            bc_name = face_mapping[rcr_id]
            vals = bc_map[bc_name]
            
            rcrt_file.write('2\n')
            rcrt_file.write(str(vals['Rp']) + '\n')
            rcrt_file.write(str(vals['C']) + '\n')
            rcrt_file.write(str(vals['Rd']) + '\n')
            rcrt_file.write('0.0 ' + str(vals['Pd']) + '\n')
            rcrt_file.write('1.0 ' + str(vals['Pd']) + '\n')
    