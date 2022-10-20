import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n
from vtk.util.numpy_support import numpy_to_vtk as n2v
import numpy as np
from .file_io import parse_mdl



class Polydata():
    
    def __init__(self):
        self.polydata = None
    
    def load_polydata(self, input_file, format = 'xml'):
        '''loads a PolyData file from <input_file> which is in <format>
        
        format: (str) One of 'xml', NULL
        '''
        
        if format == 'xml':
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(input_file)
            reader.Update()
            self.polydata = reader.GetOutput()
            

    def write_polydata(self, output_file, format = 'xml'):
        ''' writes the polydata information to <output_file> in <format>
        
        format: (str) One of 'xml', NULL
        '''
        
        if format == 'xml':
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(output_file)
            writer.SetInputData(self.polydata)
            writer.Update()
            writer.Write()

    def get_pointdata(self):
        ''' get pointdata object
        '''
        return self.polydata.GetPointData()
    
    def get_celldata(self):
        ''' get celldata object
        '''
        return self.polydata.GetCellData()

    def get_points(self):
        ''' retrieve xyz points
        '''
        return v2n(self.polydata.GetPoints().GetData())    
    
    def get_pointdata_array(self, array_name):
        ''' retrieve array data from Point Data
        '''
        return v2n(self.polydata.GetPointData().GetArray(array_name))
    
    def get_celldata_array(self, array_name):
        ''' retrieves array data from Cell Data
        '''
        return v2n(self.polydata.GetCellData().GetArray(array_name))
    
    def get_pointdata_arraynames(self):
        ''' get pointdata array names
        '''
        pointdata = self.get_pointdata()
        array_num = pointdata.GetNumberOfArrays()
        array_names = [pointdata.GetArrayName(i) for i in range(array_num)]
        return array_names
    
    def add_pointdata(self, array: np.array, array_name):
        ''' Adds a new array to point data.
        Use this to overwrite old arrays or add new arrays. To overwrite, simply provide the old array's array name.
        '''
        new_celldata = n2v(array)
        new_celldata.SetName(array_name)
        self.polydata.GetPointData().AddArray(new_celldata)
       
    def remove_pointdata_array(self, array_name):
        ''' removes an array from pointdata
        '''
        self.polydata.GetPointData().RemoveArray(array_name)

        
    def rename_pointdata_array(self, old_name, new_name):
        ''' Rename point data array
        '''
        tmp = self.polydata.GetPointData().GetArray(old_name)
        tmp.SetName(new_name)
    

class Centerlines(Polydata):
    ''' Handles centerlines
    '''
    
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
        super().__init__()

    @classmethod
    def from_file(cls, centerlines_file):
        c = cls()
        c.read_polydata(centerlines_file)
        return c
        
    def check_centerlines_data(self):
        """ Check that the centerline data contains all of the required fields.
        """
        field_names = [v for k, v in self.CellDataFields.__dict__.items() if not k.startswith('__')]
        for field in field_names:
            if not self.polydata.GetCellData().GetArray(field):
                print("Centerlines do not contain the '%s' data field." % field)
                return False

        field_names = [v for k, v in self.PointDataFields.__dict__.items() if not k.startswith('__')]
        for field in field_names:
            if not self.polydata.GetPointData().GetArray(field):
                print("Centerlines do not contain the '%s' data field." % field)
                return False
        return True

    
    def generate_centerlines(self, mdl, vtp, inlet, use_entire_tree = True, outlet_names = []):
        ''' Generates centerlines 
        
        mdl:    (str) file path to the .mdl file
        vtp:    (vtp) file path to the model .vtp file
             
        '''
    
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
        # generate centerlines for entire model
        if use_entire_tree:
            inlet_ids = [face_mappings[inlet]]
            del face_mappings[inlet]
            outlet_ids = list(face_mappings.values())
        # only generate centerlines for part of the model.
        else:
            inlet_ids = [face_mappings[inlet]]
            outlet_ids = [face_mappings[name] for name in outlet_names]

        self.polydata = sv.vmtk.centerlines(model_polydata, inlet_ids, outlet_ids, use_face_ids=True)
    