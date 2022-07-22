from src.file_io import Solver0D
import sys
if __name__ == '__main__':
    
    t = Solver0D()
    t.read_solver_file(sys.argv[1])
    
    for j in t.junctions:
        if j['junction_type'] == 'BloodVesselJunction':
            
            j['junction_type'] = 'NORMAL_JUNCTION'
            
    t.write_solver_file('temp.in')