
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.polydata import Polydata
#! Currently not working


if __name__ == '__main__':
    # Create a modeler.

    # Read model geometry.
    p = Polydata()
    # requires parasolid plugin
    p.convert_from_parasolid('data/stenosis/AS1_SU0308/Models/SU0308_prestent.xmt_txt')
    p.write_polydata('test.vtp')
    
