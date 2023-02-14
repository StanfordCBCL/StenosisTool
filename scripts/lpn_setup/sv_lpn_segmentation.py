# File: sv_lpn_segmentation.py
# File Created: Monday, 31st October 2022 7:20:33 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 13th February 2023 9:18:37 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Utilize Simvascular to construct a 0D LPN from centerlines using Martin's automated pipeline.

try:
    import sv_rom_simulation
except ImportError as e:
    print(e,': please run using simvascular --python -- this_script.py')
    exit(1)
    
from pathlib import Path
import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = "Construct a basic LPN in the workspace.")
    parser.add_argument("-mn", help = 'model_name')
    parser.add_argument("-ofc", help = 'outlet_face_names file')
    parser.add_argument("-cent", help = "centerlines file")
    parser.add_argument("-flow", help = 'Flow file')
    parser.add_argument("-ts", type = float, help = "time step")
    parser.add_argument("-n_ts", type = int, help = "number of time steps")
    parser.add_argument("-odir", help = "Outlet directory")
    
    
    
    args = parser.parse_args()

    
    ## needed to write 0d file
    params = sv_rom_simulation.Parameters()

    ## FILL PARAMS, leaving defaults
    params.model_order = 0
    params.model_name = args.mn
    params.outlet_face_names_file = args.ofc
    params.centerlines_input_file = args.cent
    
    # fluid
    params.density = 1.06
    params.viscosity = .04
    
    # material (only allows linear for now.)
    params.material_model = 'LINEAR'
    params.linear_material_ehr = 1.2e6
    params.linear_material_pressure = 0
    params.uniform_material = True
    
    # retrieve appropriate inflow
    params.inflow_input_file = args.flow
    
    # simulation params
    params.time_step = args.ts
    params.num_time_steps = args.n_ts
    
    # other
    params.output_directory = args.odir
    params.compute_centerlines = False
    params.solver_output_file = args.mn + '.in'
    
    # BC
    params.uniform_bc = False
    params.outflow_bc_type = ['rcrt.dat']
    params.outflow_bc_file = args.odir

    ## GET CENTERLINES
    centerlines = sv_rom_simulation.Centerlines()
    centerlines.read(None, params.centerlines_input_file)
    
    ## GENERATE MESH
    mesh_0d = sv_rom_simulation.Mesh()
    mesh_0d.generate(params, centerlines)