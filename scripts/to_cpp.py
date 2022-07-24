from src.solver import Solver0D
from src.misc import create_parser
import sys

if __name__ == '__main__':
    
    x = Solver0D()
    x.read_solver_file(sys.argv[1])
    x.to_cpp()
    x.write_solver_file(sys.argv[1])