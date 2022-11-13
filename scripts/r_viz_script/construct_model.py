

import numpy as np
    
import argparse
from src.polydata import Centerlines
from src.file_io import read_json
from src.misc import create_tool_parser, get_solver_path
from src.stenosis import rp_to_radius
from src.lpn import Solver0D
import numpy as np
import re
from pathlib import Path





def main(args):
    
    for solver_dir in args.solver_dirs:
        solver_dir = Path(solver_dir)
        
        # centerlines
        centerlines = Centerlines()
        centerlines.load_centerlines(solver_dir / 'model_centerlines.vtp')
        # stenosis file
        stenosis_info = read_json(solver_dir / 'stenosis_vessels.dat')
        # solver file
        solver_file = get_solver_path(solver_dir)
        solver0d = Solver0D()
        solver0d.read_solver_file(solver_file)        
        
        # retrieve branch ids and create stenosis array
        branch_ids = centerlines.get_pointdata(Centerlines.PointDataFields.BRANCHID)
        stenosis_nodes = np.ones_like(branch_ids) * -1

        mu = solver0d.simulation_params['viscosity']
        # for each branch that we are modifying
        for idx, vessels in enumerate(stenosis_info['all_changed_vessels']):
            # find branch id of vessels
            v = solver0d.get_vessel(vessels[0])
            branch_id = re.search("branch([0-9]+)_.*", v['vessel_name']).group(1)
            stenosis_nodes = np.where(branch_ids == int(branch_id), idx, stenosis_nodes )
            
            #! THIS IS STILL MISSING LOGIC FOR CASES WHERE YOU MODIFY THE MIDDLE SEGMENT ONLY/A PARTIAL OF THE WHOLE BRANCH
            r = rp_to_radius(vessels, mu)
            
            












def get_ref_frame(vec):
    """
    Generate a reference frame given an axial vector vec
    """
    v1 = (vec)/np.linalg.norm(vec)
    v2 = np.cross(v1,np.array([0,0,1]))
    if(np.linalg.norm(v2) > 1.0e-4):
        v2 = np.cross(v1,np.array([0,1,0]))
    v2 = v2/np.linalg.norm(v2)
    v3 = np.cross(v1,v2)
    v3 /= np.linalg.norm(v3)
    return np.concatenate((v1.reshape(1,-1),v2.reshape(1,-1),v3.reshape(1,-1)),axis=0)

def get_points(target, origin):
    """
    Create all ponts for the cylindrical geometry
    """
    # Get segment length
    segment_lengh = np.linalg.norm(target-origin)
    # Create reference frams
    # ref[0] is the axial direction 
    # ref[1] and ref[2] are the two orthogonal directions
    ref = get_ref_frame(target-origin)
    # Create points 
    points = np.zeros((n_segments*n_theta,3))
    for loopA in range(n_segments):
        curr_center = origin + loopA*(segment_lengh/n_segments)*ref[0]
        for loopB in range(n_theta):
            points[loopA*n_theta + loopB] = curr_center + ref[1]*radius*np.cos(2.0*np.pi*loopB/n_theta) + ref[2]*radius*np.sin(2.0*np.pi*loopB/n_theta)
    return points

# MAIN CODE
if __name__ == "__main__":
    
    # parse for centerlines
    parser = create_tool_parser('Add stenosis locations')
    args = parser.parse_args()
    main(args)
    
    
    
    
    
    
    
    
    

    # Set parameters
    # Number of axial segments
    n_segments = 10
    # Number of circumferential points
    n_theta = 10
    # Cylinder radius
    radius = 10.0
    # Start point for the axial vector
    origin = np.array([0,0,0])
    # End point for the axial vector
    target = np.array([10.0,10.0,10.0])

    # Create VTK file header
    f = open("test_file.vtk", "w")
    f.write('# vtk DataFile Version 2.0\n')
    f.write('Really cool data\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')

    # Generate node coordinates
    nodeCoords = get_points()

    # Write points to file
    f.write('POINTS %d float \n' % (n_theta*n_segments))
    for loopA in range(len(nodeCoords)):
        f.write('%e %e %e \n' % (nodeCoords[loopA,0],nodeCoords[loopA,1],nodeCoords[loopA,2]))

    # Separate
    f.write('\n')

    # Write triangle stript to file
    f.write('TRIANGLE_STRIPS %d %d\n' % (n_segments-1,(n_segments-1)*(2*n_theta+3)))
    for loopA in range(n_segments-1):
        f.write('%d ' % (2*n_theta+2))
    for loopB in range(n_theta):
        f.write('%d %d ' % (loopA*n_theta+loopB,(loopA+1)*n_theta+loopB))
    f.write('%d %d\n' % (loopA*n_theta,(loopA+1)*n_theta))

    f.close()