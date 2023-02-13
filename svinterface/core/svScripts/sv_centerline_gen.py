try:
    import sv
except ImportError as e:
    print(e, ': please run using simvascular --python -- this_script.py')
    exit(1)
    
import xml.etree.ElementTree as ET
import sys
import vtk
    
    
def parse_mdl(mdl_file: str, reverse = False):
    """parses mdl file for cap name mapping to faceid

    Args:
        mdl_file (str): path to mdl file
        reverse (bool, optional): reverses mapping so faceid maps to cap name. Defaults to False.

    Returns:
        dict: a dict describing map of cap name to faceid
    """
    
    # since .mdl has a line like "<format version="1.0" />" which fails for the standard XMLparser, rather than creating a custom parser, just remove that line after reading the file in and parse as a list of strings
    mdl_file = str(mdl_file)
    with open(mdl_file, 'r') as mdl:
        lines = mdl.readlines()
        if 'format' in lines[1]:
            lines[1] = ''
    
    root = ET.fromstringlist(lines)
    faces = root[0][0].find('faces')

    # save to a dict
    face_mappings = {}
    for face in faces:
        if face.attrib['type'] == 'cap':
            if reverse:
                face_mappings[int(face.attrib['id'])] = face.attrib['name']
            else:
                face_mappings[face.attrib['name']] = int(face.attrib['id'])
    return face_mappings


def generate_centerlines( vtp, inlet_face_ids, outlet_face_ids):
    ''' Generates centerlines
    
    vtp:    (vtp) file path to the model .vtp file
    outlet_face_ids: face ids to generate centerlines to.
    '''
    
    # Create a modeler.
    kernel = sv.modeling.Kernel.POLYDATA
    modeler = sv.modeling.Modeler(kernel)

    # Read model geometry.
    print('Generating centerlines from:', vtp)
    model = modeler.read(vtp)
    model_polydata = model.get_polydata()

    return sv.vmtk.centerlines(model_polydata, inlet_face_ids, outlet_face_ids, use_face_ids=True)


if __name__ == '__main__':
    
    # parse the input from subprocess
    vtp = sys.argv[1]
    inlet_ids = [int(sys.argv[2])]
    outlet_ids = [int(idx) for idx in sys.argv[3].split("|")]
    outfile = sys.argv[4]
    
    c = generate_centerlines(vtp,inlet_face_ids=inlet_ids, outlet_face_ids=outlet_ids)
    
    # write centerlines
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outfile)
    writer.SetInputData(c)
    writer.Update()
    writer.Write()
