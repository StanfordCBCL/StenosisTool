# File: add_solver_rcrs.py
# File Created: Thursday, 28th July 2022 11:21:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 5th August 2022 1:36:51 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: If an rcrt and a solver file exist in a dir, maps the solver rcrt file into the solver file



from src.misc import create_tool_parser, get_solver_name
from src.bc import BoundaryConditions
from src.solver import Solver0D
import os


def add_rcrs(solver: Solver0D, bcs: BoundaryConditions):
    ''' Adds appropriate rcrs to the dummy solver 
    '''
    
    bc_inv_map = bcs.get_bc_map()
    
    for bc in solver.bc:
        if bc['bc_type'] == 'RCR':
            bc_values = bc_inv_map[solver.bc_map[bc['bc_name']]]
            bc['bc_values']['Rp'] = bc_values['Rp']
            bc['bc_values']['Rd'] = bc_values['Rd']
            bc['bc_values']['C'] = bc_values['C']
            bc['bc_values']['Pd'] = bc_values['Pd']
    


if __name__ == '__main__':
    
    tool_parser = create_tool_parser('Adds existing rcr file to a solver file')
    
    args = tool_parser.parse_args()
    
    for sdir in args.solver_dirs:
        solver_file = os.path.join(sdir, get_solver_name(sdir))
    
        solver = Solver0D()
        solver.read_solver_file(solver_file)
        
        bc = BoundaryConditions()
        bc.read_rcrt_file(rcrt_file=os.path.join(sdir, 'rcrt.dat' ))
        
        add_rcrs(solver, bc)