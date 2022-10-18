# File: add_solver_rcrs.py
# File Created: Thursday, 28th July 2022 11:21:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 18th October 2022 12:03:00 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: If an rcrt and a solver file exist in a dir, maps the solver rcrt file into the solver file


from src.parser import ToolParser
from src.data_org import LPNDir
from src.bc import BoundaryConditions
from src.lpn import LPN


def add_rcrs(lpn: LPN, bcs: BoundaryConditions):
    ''' Adds appropriate rcrs to the dummy solver 
    '''
    
    bc_inv_map = bcs.get_bc_map()
    
    for bc in lpn.bc:
        if bc['bc_type'] == 'RCR':
            bc_values = bc_inv_map[lpn.bc_map[bc['bc_name']]]
            bc['bc_values']['Rp'] = bc_values['Rp']
            bc['bc_values']['Rd'] = bc_values['Rd']
            bc['bc_values']['C'] = bc_values['C']
            bc['bc_values']['Pd'] = bc_values['Pd']
    


if __name__ == '__main__':
    
    tool_parser = ToolParser('Adds existing rcr file to a lpn file')
    
    args = tool_parser.parse_args()
    
    for lpn_dir in args.lpn_dirs:
        
        # writes existing rcr into lpn file
        lpn_dir = LPNDir(lpn_dir)
    
        lpn = LPN.from_file(str(lpn_dir.lpn_path))
        
        bc = BoundaryConditions()
        bc.read_rcrt_file(rcrt_file=str(lpn_dir.rcrt_file))
        
        add_rcrs(lpn, bc)
        
        lpn.write_lpn_file(str(lpn_dir.lpn_path))