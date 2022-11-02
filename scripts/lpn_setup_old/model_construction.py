# File: jc_model_construction.py
# File Created: Monday, 18th July 2022 3:40:19 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 27th October 2022 5:03:06 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Constructs a new directory from an base one where everything is identical except the solver file is processed to includes junction coefficients


from src.lpn import LPN
from src.data_org import LPNDir
from src.parser import ToolParser
import shutil
from pathlib import Path

def convert_to_jc(lpn: LPN):
    
    for junc in lpn.junctions:
        if junc['junction_type'] != 'internal_junction':
            
            inlet = junc['inlet_vessels']
            outlets = junc['outlet_vessels']
            S0 = junc['areas'][0]
            s_outlets = junc['areas'][1:]
            
            assert len(outlets) == len(s_outlets), 'Number of areas does not match number of outlets'
            
            outlet_areas = list(zip(outlets, s_outlets))
            
            density = lpn.simulation_params['density']
            for out, S1 in outlet_areas:
                out_vess = lpn.get_vessel(out)
                
                # update with stenosis coefficient of junction
                out_vess['junction_coefficient'] = 1.52 * density * ( (S0/S1 - 1) **2) / (2 * S0**2)
                out_vess['zero_d_element_values']['stenosis_coefficient'] += out_vess['junction_coefficient']

def main(args):
    print(f'Constructing LPN dir for {args.out_dir}...',end = '\t', flush = True)
    

    
    out_dir = Path(args.out_dir)
    # construct a copy of the base directory
    if out_dir.exists() and not args.force:
        print(str(out_dir) + ' already exists. Skipping')
        return
    
    # copy all base files
    out_dir.mkdir(exist_ok = True)
    for file in Path(args.in_dir).iterdir():
        shutil.copy(str(file), str(out_dir))
    print('Done')
    
    # copy over centerlines
    
    if args.jc:
        print('Converting to JC LPN...', end = '\t', flush = True)
        lpn_dir = LPNDir(args.out_dir)
        lpn_jc = LPN.from_file(str(lpn_dir.lpn_path))
        convert_to_jc(lpn_jc)
        new_jc_lpn = lpn_dir.lpn_path.with_suffix('.jc.in')
        lpn_dir.lpn_path.rename(new_jc_lpn)
        lpn_jc.write_lpn_file(str(new_jc_lpn))
        print('Done')
    

        

if __name__ == '__main__':
    
    tool = Parser(desc = 'Set up a model')
    
    tool.parser.add_argument('-model', dest = 'lpn_di', help = 'Dirname of base files to copy' )
    tool.parser.add_argument('-i', dest = 'in_dir', help = 'Dirname of base files to copy' )
    tool.parser.add_argument('-o', dest = 'out_dir', default = 'base_lpn_dir', help = 'Dirname of output LPN dir')
    tool.parser.add_argument('-jc', dest = 'jc', default = False, action = 'store_true', help = 'Flag to convert LPN file to a JC LPN file.')
    tool.parser.add_argument('-f', dest = 'force', default = False, action = 'store_true', help = 'Force override files in existing dir')
    args = tool.parse_args()

    main(args)
    
    