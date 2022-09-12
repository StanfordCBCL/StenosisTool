# File: sv_dev_segmentation_gen.py
# File Created: Thursday, 28th July 2022 3:33:40 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Friday, 5th August 2022 1:28:14 am
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Using Simvascular Automated Pipeline & Development given parameters, construct a 0D model from a particular 3D geometry



try:
    import sv_rom_simulation
except ImportError as e:
    print(e,': please run using simvascular --python -- this_script.py')
    exit(1)

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.misc import *
from src.flow import Inflow0D
from src.file_io import parse_mdl
from src.data_org import DataPath

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
    def __init__(self, model_name, inlet, outlets):
        self.model_name = model_name
        self.model_order = 0 
        self.inlet = inlet
        self.outlets = outlets
        self.units = 'mm'
        
class FileParams():
    def __init__(self, centerlines_file, solver_dir, inflow_file, out_name):
        self.solver_dir = solver_dir
        self.centerlines_input_file = centerlines_file
        self.inflow_file = inflow_file
        self.out_name = out_name
        self.inverse_flow = False
        self.smooth = True
        
class SimulationParams():
    def __init__(self):
        self.num_cycles = 6
        self.num_ts_per_cycle = 1000
    

def write_0d_file(files: FileParams, model: ModelParams, mesh: MeshParams, fluids: FluidParams, material: MaterialParams, sim: SimulationParams, dummy_bc = True):
    ''' Writes a 0D file '''
    
    ## needed to write 0d file
    params = sv_rom_simulation.Parameters()

    ## params
    
    params.model_order = model.model_order
    params.model_name = model.model_name
    
    with open(os.path.join(files.solver_dir, 'inlet_face_names.dat'), 'w') as infile:
        infile.write(model.inlet + '\n')
    with open(os.path.join(files.solver_dir, 'outlet_face_names.dat'), 'w') as outfile:
        for outlet in model.outlets:
            outfile.write(outlet + '\n')

    params.outlet_face_names_file = os.path.join(files.solver_dir, 'outlet_face_names.dat')
    params.centerlines_input_file = files.centerlines_input_file
    

    # mesh params
    params.element_size = mesh.element_size
    params.seg_size_adaptive = mesh.seg_size_adaptive
    params.seg_size = mesh.seg_size
    params.seg_min_size = mesh.seg_min_size
    
    # fluid
    params.density = fluids.density
    params.viscosity = fluids.viscocity
    
    # material
    params.material_model = 'LINEAR'
    params.linear_material_ehr = material.linear_ehr
    params.linear_material_pressure = material.linear_pressure
    params.uniform_material = True
    

    # change Units
    params.set_units(model.units)
    
    
    inflow = Inflow0D.from_file(files.inflow_file, inverse = files.inverse_flow, smooth = files.smooth)
    smooth_flowfile = os.path.join(files.solver_dir, 'inflow.flow')
    inflow.write_flow(smooth_flowfile)
    
    # simulation params
    params.time_step = inflow.tc / sim.num_ts_per_cycle
    params.num_time_steps = sim.num_ts_per_cycle * (sim.num_cycles - 1)
    
    # other
    params.output_directory = files.solver_dir
    params.compute_centerlines = False
    params.solver_output_file = files.out_name
    
    
    mesh_0d = sv_rom_simulation.Mesh()
    
    # put directory in outflowbcfile
    # put bc type in outflow_bc_Type as list
    params.uniform_bc = False   
    params.outflow_bc_type = ['rcrt.dat']
    
    rcr_file = os.path.join(files.solver_dir, 'rcrt.dat')
    if dummy_bc:
        # just write an empty file
        if not os.path.exists(rcr_file):
            with open(rcr_file, 'w') as rfile:
                rfile.write('2\n') # dummy rcrt
                for outlet in model.outlets:
                    rfile.write('2\n')
                    rfile.write(outlet + '\n')
                    rfile.write('0\n')
                    rfile.write('0\n')
                    rfile.write('0\n')
                    rfile.write('0 0\n')
                    rfile.write('0 0\n')
    
    
    params.outflow_bc_file = files.solver_dir
    params.inflow_input_file = smooth_flowfile


    # centerlines for mesh
    centerlines = sv_rom_simulation.Centerlines()
    centerlines.read(params, params.centerlines_input_file)
    
    mesh_0d.generate(params, centerlines)
    
    return



########
# Main #
########

def main(args):
    
    org = DataPath(args.root)
    if args.models:
        model_list = args.models
    else:
        model_list = list(org.model_names)
        
    for model_name in model_list:
        print('Retrieving segmentation for {}...'.format(model_name), end = '\t', flush = True)
        model = org.find_model(model_name)
        
        
        #! Can implement a config file here to change params
        # model
        face_mappings = parse_mdl(model.info['files']['mdl_file'])
        del face_mappings[model.info['model']['inlet']]
        outlets = face_mappings.keys()
        mod_params = ModelParams(model_name=model.info['metadata']['name'],
                                 inlet=model.info['model']['inlet'],
                                 outlets=outlets)
    
        mod_params.units = model.info['model']['units']
        # mesh
        mesh_params = MeshParams()

        # simulation
        sim_params = SimulationParams()
        
        # fluid
        fluid_params = FluidParams()

        # material
        material_params = MaterialParams()
        
        # files
        if model.info['files']['rom_inflow']:
            inflow = model.info['files']['rom_inflow']
            inverse = False
        elif model.info['files']['inflow']:
            inflow = model.info['files']['inflow']
            inverse = True
        else:
            raise RuntimeError('Inflow file does not exist')

        file_params = FileParams(centerlines_file=model.model_centerlines,
                                 solver_dir=model.base_solver_dir,
                                 inflow_file=inflow,
                                 out_name=os.path.basename(model.base_model_solver))
        file_params.inverse_flow = inverse
        file_params.smooth = True

        write_0d_file(files=file_params,
                      model=mod_params,
                      mesh=mesh_params,
                      fluids=fluid_params,
                      material=material_params,
                      sim=sim_params)
        

        print('Done')
           
             
        

if __name__ == '__main__':
    
    dev_parser = create_dev_parser(desc= 'Constructs a 0d input segmentation of the vasculature without BC for a 0d solver for pulmonary vasculature')
    
    args = dev_parser.parse_args()
    main(args)
        