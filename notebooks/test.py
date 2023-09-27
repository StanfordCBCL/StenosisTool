from svinterface.core.polydata import Polydata, LegacyVTK
import vtk
from vtk.util.numpy_support import numpy_to_vtk as n2v
from vtk.util.numpy_support import vtk_to_numpy as v2n

x = Polydata.load_polydata('data/diseased/AS1_SU0308_stent/Segmentations/RPA_stented_stented.vtp')

data = vtk.vtkIdTypeArray()
x.polydata.GetLines().ExportLegacyFormat(data)
print(v2n(data))
y = LegacyVTK()

y.add_polydata(x.get_points(), v2n(data))

y.write_vtk("images/test.vtk")
