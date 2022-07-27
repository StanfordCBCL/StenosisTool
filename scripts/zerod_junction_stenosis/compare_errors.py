from src.solver import Solver0D
from src.solver_results import SolverResults
from src.centerlines import Centerlines
from src.misc import get_basename
import argparse


def get_caps():
    # TODO: Use the RCR map in the solver file, iterate through the tree, and obtain the branch numbers of the caps as well as their vessel ID. Then take the vtp file for the 3D and identify the last point in the branch, Return a list of their two locations, in pairs
    
    pass

def compute_errors():
    #TODO: Retrieving pressure/flow values at particular points in particular arrays to compare values at caps
    pass

def main(args):
    # files
    solver_file = args.solver_file
    solver_result_csv = get_basename(solver_file) + '_branch_results.csv'
    three_d = args.three_d
    
    solver0D = Solver0D()
    solver0D.read_solver_file(solver_file)
    
    results = SolverResults.load_from_csv()
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest = 'solver_file', help = '0D solver file containing vasculature to compute error for')
    parser.add_argument('-t', dest = 'three_d', help= '1D mapped representation of the 3D model')

    args = parser.parse_args()
    
    main(args)