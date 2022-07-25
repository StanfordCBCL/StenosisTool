

import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n
from vtk.util.numpy_support import numpy_to_vtk as n2v
import numpy as np
from .file_io import parse_mdl


class Centerlines():
    ''' Handles centerlines '''
    class CellDataFields(object):
        pass

    class PointDataFields(object):
        """ This class defines the standard used field point data field names.
        """
        AREA = "CenterlineSectionArea"
        CENTID = "CenterlineId"
        PATH = "Path"
        BRANCHID = "BranchId"
        BIFURCATIONID = "BifurcationId"
        NODEID = "GlobalNodeId"
        NORMAL = "CenterlineSectionNormal"
        
    def __init__(self):
        self.centerlines = None
        
    def load_centerlines(self, input_file):
        '''loads a centerlines file'''
        # Create a modeler.
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(input_file)
        reader.Update()
        self.centerlines = reader.GetOutput()
    
    def write_centerlines(self, output_file):
        # write centerline to file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_file)
        writer.SetInputData(self.centerlines)
        writer.Update()
        writer.Write()
        
    def check_centerlines_data(self):
        """ Check that the centerline data contains all of the required fields.
        """
        field_names = [v for k, v in self.CellDataFields.__dict__.items() if not k.startswith('__')]
        for field in field_names:
            if not self.centerlines.GetCellData().GetArray(field):
                print("Centerlines do not contain the '%s' data field." % field)
                return False

        field_names = [v for k, v in self.PointDataFields.__dict__.items() if not k.startswith('__')]
        for field in field_names:
            if not self.centerlines.GetPointData().GetArray(field):
                print("Centerlines do not contain the '%s' data field." % field)
                return False

        return True
    
    def get_pointdata(self, array_name):
        ''' retrieve array data from Point Data'''
        return v2n(self.centerlines.GetPointData().GetArray(array_name))
    
    def get_celldata(self, array_name):
        ''' retrieves array data from Cell Data'''
        return v2n(self.centerlines.GetCellData().GetArray(array_name))

    def get_pointdata_arraynames(self):
        array_num = self.centerlines.GetPointData().GetNumberOfArrays()
        array_names = [self.centerlines.GetPointData().GetArrayName(i) for i in range(array_num)]
        return array_names
    
    def add_pointdata(self, array: np.array, array_name):
        ''' Use this to overwrite old arrays or add new arrays. To overwrite, simply provide the old array's array name'''
        new_celldata = n2v(array)
        new_celldata.SetName(array_name)
        self.centerlines.GetPointData().AddArray(new_celldata)
    
    def rename_pointdata(self, old_name, new_name):
        tmp = self.centerlines.GetPointData().GetArray(old_name)
        tmp.SetName(new_name)
    
    def remove_pointdata(self, array_name):
        self.centerlines.GetPointData().RemoveArray(array_name)
    
    def generate_centerlines(self, mdl, vtp, inlet, use_entire_tree = True, outlet_names = []):
        ''' Generates centerlines '''
    
        try:
            import sv
        except ImportError as e:
            print(e, ': please run using simvascular --python -- this_script.py')
            exit(1)
        
        # Create a modeler.
        kernel = sv.modeling.Kernel.POLYDATA
        modeler = sv.modeling.Modeler(kernel)

        # Read model geometry.
        model = modeler.read(vtp)
        
        # parse mdl for mapping
        face_mappings = parse_mdl(mdl)
        model_polydata = model.get_polydata()

        ## Calculate centelines. 
        # Use face IDs.
        if use_entire_tree:
            inlet_ids = [face_mappings[inlet]]
            del face_mappings[inlet]
            outlet_ids = list(face_mappings.values())
        else:
            inlet_ids = [face_mappings[inlet]]
            outlet_ids = [face_mappings[name] for name in outlet_names]

        self.centerlines = sv.vmtk.centerlines(model_polydata, inlet_ids, outlet_ids, use_face_ids=True)
    