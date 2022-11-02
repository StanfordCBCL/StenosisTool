# File: sv_lpn_segmentation.py
# File Created: Monday, 31st October 2022 7:20:33 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Tuesday, 1st November 2022 4:08:28 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Construct a basic LPN in the workspace.

try:
    import sv_rom_simulation
except ImportError as e:
    print(e,': please run using simvascular --python -- this_script.py')
    exit(1)
    
import sys
import os
from pathlib import Path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


from sgt.utils.parser import ToolParser
import sgt.utils.io as io
from sgt.core.manager import LPNConstructionManager
from sgt.core.flow import Inflow
from sgt.core.lpn import LPN
from sgt.core.bc import BoundaryConditions

##########
# Params #
##########

class MaterialParams():
    def __init__(self) -> None:
        self.linear_ehr = 1.2e6 # dyne/cm^2
        self.linear_pressure = 0

class FluidParams():
    def __init__(self):
        self.density = 1.06
        self.viscocity = 0.04

class MeshParams():
    def __init__(self):
        self.element_size = .1
        self.seg_size_adaptive = False
        self.seg_size = 99999
        self.seg_min_size = 1
        
class ModelParams():
    def __init__(self, model_name, inlet, outlets, units = 'cms'):
        self.model_name = model_name
        self.model_order = 0 
        self.inlet = inlet
        self.outlets = outlets
        self.units = units
        
class SimulationParams():
    def __init__(self):
        self.num_cycles = 6
        self.num_ts_per_cycle = 1000

########
# Code #
########

def segment_lpn(M: LPNConstructionManager, model: ModelParams, mesh: MeshParams, fluids: FluidParams, material: MaterialParams, sim: SimulationParams):
    ''' Writes a 0D LPN file 
    '''
    
    ## needed to write 0d file
    params = sv_rom_simulation.Parameters()

    ## FILL PARAMS
    params.model_order = model.model_order
    params.model_name = model.model_name
    
    # write inlet/outlet to a dat file
    with (M.lpn_files / 'inlet_face_names.dat').open('w') as infile:
        infile.write(model.inlet + '\n')
    with (M.lpn_files / 'outlet_face_names.dat').open('w') as outfile:
        for outlet in model.outlets:
            outfile.write(outlet + '\n')

    params.outlet_face_names_file = str(M.lpn_files / 'outlet_face_names.dat')
    params.centerlines_input_file = str(M.centerlines)
    

    # mesh params
    params.element_size = mesh.element_size
    params.seg_size_adaptive = mesh.seg_size_adaptive
    params.seg_size = mesh.seg_size
    params.seg_min_size = mesh.seg_min_size
    
    # fluid
    params.density = fluids.density
    params.viscosity = fluids.viscocity
    
    # material (only allows linear for now.)
    params.material_model = 'LINEAR'
    params.linear_material_ehr = material.linear_ehr
    params.linear_material_pressure = material.linear_pressure
    params.uniform_material = True
    
    # change Units
    params.set_units(model.units)
    
    # retrieve appropriate inflow
    inflow = Inflow.from_file(str(M.flow), smooth = True)
    smooth_flowfile = str(M.lpn_files / 'smooth_inflow.flow')
    inflow.write_flow(smooth_flowfile)
    params.inflow_input_file = smooth_flowfile
    
    # simulation params
    params.time_step = inflow.tc / sim.num_ts_per_cycle
    params.num_time_steps = sim.num_ts_per_cycle * (sim.num_cycles - 1)
    
    # other
    params.output_directory = str(M.lpn_files)
    params.compute_centerlines = False
    params.solver_output_file = str(M.lpn.name)
    
    # BC
    params.uniform_bc = False
    params.outflow_bc_type = ['rcrt.dat']
    
    # tune mode
    if M.tune:
        rcrt_path = Path(M.lpn_files.name) / 'rcrt.dat'
        # add to config
        M.config_add(['paths','rcrt_file'], str(rcrt_path))
        M.write_config()
        
        # write an empty if BC does not exist
        bc = BoundaryConditions()
        for outlet in model.outlets:
            bc.add_rcr(outlet, 0, 0, 0, 0)
        bc.write_rcrt_file(str(M.lpn_files))
        params.outflow_bc_file = str(M.lpn_files)
    else:
        if M.rcrt is None:
            raise FileNotFoundError("rcrt file not found, but BC tuning was not requested.")
        params.outflow_bc_file = str(M.rcrt.parent)


    ## GET CENTERLINES
    centerlines = sv_rom_simulation.Centerlines()
    centerlines.read(None, params.centerlines_input_file)
    
    ## GENERATE MESH
    mesh_0d = sv_rom_simulation.Mesh()
    mesh_0d.generate(params, centerlines)
    
    return


if __name__ == '__main__':
    
    parser = ToolParser(desc = "Construct a basic LPN in the workspace.")
    args = parser.parse_args()
    
    M = LPNConstructionManager(args.config)
        
    ## Fill Params (Use Defaults)
    # model
    face_ids = io.parse_mdl(M.mdl)
    del face_ids[M.inlet]
    outlets = list(face_ids.keys())
    mod_params = ModelParams(model_name=M.model_name,
                                inlet=M.inlet,
                                outlets=outlets,
                                units = M.units)

    # mesh
    mesh_params = MeshParams()

    # simulation
    sim_params = SimulationParams()
    
    # fluid
    fluid_params = FluidParams()

    # material
    material_params = MaterialParams()

    print("Segmenting " + M.model_name + "...", end = '\t', flush = True)
    # write the LPN file
    segment_lpn(M=M,
                model=mod_params,
                mesh=mesh_params,
                fluids=fluid_params,
                material=material_params,
                sim=sim_params)
    
    # read LPN
    lpn = LPN.from_file(str(M.lpn))
    
    # construct bc map
    lpn.write_bc_map(outlets)
    lpn.to_cpp()
    
    lpn.write_lpn_file(M.lpn)
    
    
    
    print('Done')