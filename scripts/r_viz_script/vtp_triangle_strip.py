import numpy as np

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

def get_points():
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

  # Set parameters
  # Number of axial segments
  n_segments = 2
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