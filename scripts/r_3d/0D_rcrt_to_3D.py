# File: 0D_rcrt_to_3D.py
# File Created: Monday, 31st October 2022 8:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 23rd January 2023 7:45:09 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Converts a 0D rcrt file to 3D rcrt file.


import re
from sgt.utils.parser import Parser
from sgt.core.bc import BoundaryConditions
from pathlib import Path

def parse_svpre(svpre: Path):
    ''' parses an svpre file to extract face mappings
    '''
    id_to_name = {}
    with svpre.open('r') as sfile:
        for line in sfile:
            match= re.search("set_surface_id_vtp mesh-complete/mesh-surfaces/(.*).vtp ([0-9]*)", line)
            if match:
                id_to_name[int(match.group(2))] = match.group(1)
                
    return id_to_name
    
def parse_inp(inp: Path):
    ''' parses an inp file to extract the order of RCRTs in the rcrt file
    '''
    with inp.open('r') as ifile:
        f = ''.join(ifile.readlines())
        rcr_num = int(re.search("^Number of RCR Surfaces: ([0-9]*)$", f, re.MULTILINE).group(1))
        
        rcr_list = (re.search("^List of RCR Surfaces: (.*)$", f, re.MULTILINE).group(1).split(' '))
        rcr_list = [int(x) for x in rcr_list]
    
    return rcr_num, rcr_list


if __name__ == '__main__':
    
    parser = Parser(desc = "Converts a 3D rcrt file to 1D rcrt file")
    parser.parser.add_argument("-rcrt", dest = "rcrt_file", help = "0D rcrt file")
    parser.parser.add_argument("-o", dest = "out_dir", help = "output dir to write rcrt.dat")
    parser.parser.add_argument("-svpre", dest = "svpre_file", help = "svpre file for 3D simulation")
    parser.parser.add_argument("-inp", dest = "inp_file", help = "inp file used in 3D simulation")
    parser.parser.add_argument("-mm", action = 'store_true', default = False, help = 'units of original 3D model')

    args = parser.parse_args()
    
    bc = BoundaryConditions()
    
    # read the 0D rcrt
    bc.read_rcrt_file(args.rcrt_file, three_d=False)
    if args.mm:
        bc.cm_to_mm()
    bc_map = bc.get_bc_map()
    
    # parse files
    rcrt_num, rcrt_list = parse_inp(Path(args.inp_file))
    id_to_name = parse_svpre(Path(args.svpre_file))
    
    assert rcrt_num == len(bc.bc_list), "Number of RCR Surfaces doesn't match number of BC's listed. "
    # order bc_list
    bc.bc_list = []
    for id in rcrt_list:
        bc.bc_list.append(bc_map[id_to_name[id]])
    
    # write file out
    bc.write_rcrt_file(dirpath = args.out_dir, three_d=True)
    