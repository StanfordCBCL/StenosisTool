try:
    import sv
except ImportError as e:
    print(e + ': use simvascular --python -- this_script.py')
    exit(1)
  
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.centerlines import Polydata

if __name__ == '__main__':
    # Create a modeler.

    kernel = sv.modeling.Kernel.PARASOLID
    print(kernel)
    modeler = sv.modeling.Modeler(kernel)

    # Read model geometry.
    p = Polydata()
    model = modeler.read('data/stenosis/AS1_SU0308/Models/SU0308_prestent.xmt_txt')
    p.polydata = model.get_polydata()
    p.write_polydata('test.xml')
    
