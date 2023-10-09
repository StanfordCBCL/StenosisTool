import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n
from vtk.util.numpy_support import numpy_to_vtk as n2v
import numpy as np
from pathlib import Path
from svinterface.utils.io import parse_mdl
import subprocess

#from .file_io import parse_mdl

class LegacyVTK():
    """Legacy VTK file to construct a basic 3D geometry using points and triangle strips
    """
    def __init__(self):
        self.points = None
        self.triangle_strips = None
        
    def write_vtk(self, filename, desc = "Really cool data"):
        ''' writes vtk file to filename'''
        with open(filename, "w") as vtkfile:
            # header
            vtkfile.write("# vtk DataFile Version 2.0\n")    
            vtkfile.write(desc + "\n")
            vtkfile.write('ASCII\n')
            
            # write polydata
            vtkfile.write("DATASET POLYDATA\n")
            
            # write points
            vtkfile.write("POINTS {} float\n".format(len(self.points)))
            for point in self.points:
                vtkfile.write("%e %e %e\n" % (point[0], point[1], point[2]))
            
            vtkfile.write("\n")
            # write cells
            vtkfile.write("TRIANGLE_STRIPS {} {}\n".format(len(self.triangle_strips), self.triangle_strips.shape[0] * self.triangle_strips.shape[1] + self.triangle_strips.shape[0] ))
            for cell in self.triangle_strips:
                vtkfile.write("{} ".format(len(cell)))
                for val in cell[:-1]:
                    vtkfile.write("{} ".format(val))
                vtkfile.write("{}\n".format(cell[-1]))
                
    def add_polydata(self, points, cells):
        if self.points is None:
            self.points = points
        else:
            # concat
            self.points = np.concatenate((self.points, points))
        
        if self.triangle_strips is None:
            self.triangle_strips = cells
        else:
            # concat
            self.triangle_strips = np.concatenate((self.triangle_strips, cells))
        

class Polydata():
    """Base Polydata Class
    """
    def __init__(self, polydata = None):
        self.polydata = polydata
    
    def convert_from_parasolid(self, parasolid_file):
        '''loads a Parasolid xmt_txt file which is converted into polydata.
        Only accessible throught a private Parasolid plugin.
        '''
        #! Currently not working as intended. Please use Simvascular's GUI to perform this operation.
        raise NotImplementedError("Currently not working as intended. Please use Simvascular's GUI to perform this operation. ")
        
        try:
            import sv
        except ImportError as e:
            print(e + ': use simvascular --python -- this_script.py')
            exit(1)
        
        try:
            kernel = sv.modeling.Kernel.PARASOLID
            modeler = sv.modeling.Modeler(kernel)

            # Read model geometry.
            model = modeler.read(parasolid_file)
            self.polydata = model.get_polydata()
        except Exception as e:
            print(e, ': requires parasolid plugin, which is a private plugin.')
    
    @classmethod
    def load_polydata(cls, input_file, format = 'xml'):
        '''loads a PolyData file from <input_file> which is in <format>
        
        format: (str) One of 'xml', NULL
        '''
        
        if format == 'xml':
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(input_file)
            reader.Update()
            return cls(reader.GetOutput())
        
    def read_polydata(self, input_file):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(input_file)
        reader.Update()
        self.polydata = reader.GetOutput()
    
    @classmethod
    def create_new(cls):
        """Creates a new polydata
        """
        return cls(vtk.vtkPolyData())

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
        
    def __init__(self, centerlines = None):
        super().__init__(centerlines)
        self.centerlines = self.polydata

    @classmethod
    def load_centerlines(cls, centerlines_file):
        return cls.load_polydata(centerlines_file)
        
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

    
    def generate_centerlines(self, mdl, vtp, inlet, outfile):
        ''' Generates centerlines by creating a subprocess and calling sv.
        
        mdl:    (str) file path to the .mdl file
        vtp:    (vtp) file path to the model .vtp file
        
             
        '''
        
                
        # parse mdl for mapping
        face_mappings = parse_mdl(mdl)

        # get face ids.
        inlet_ids = [face_mappings[inlet]]
        del face_mappings[inlet]
        ordered_outlets = sorted(list(face_mappings.keys()))
        outlet_ids = [face_mappings[key] for key in ordered_outlets] #! there is an unavoidable ordering bug, so this fixes it.
        
        
        try:
            x = subprocess.run(["simvascular", "--python", "--",  str(Path(__file__).parent / "svScripts" / "sv_centerline_gen.py"), vtp, str(inlet_ids[0]), '|'.join([str(idx) for idx in outlet_ids]), outfile])
            if x.returncode!= 0:
                raise Exception("Centerlines could not be generated.")
        except FileNotFoundError:
            print("Simvascular is likely not installed on your machine.")
        
        self.polydata = self.read_polydata(outfile)
        
        return ordered_outlets
        
        
        
        
        