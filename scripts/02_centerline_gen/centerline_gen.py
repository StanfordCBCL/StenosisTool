# File: centerline_gen.py
# File Created: Sunday, 19th February 2023 11:33:21 am
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 13th July 2023 1:53:55 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Centerline Generation for prestent model / main model provided in config. Automatically tracks and adds to workspace

from svinterface.core.polydata import Centerlines
from svinterface.manager import Manager
import argparse
from pathlib import Path

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Centerline Generation for prestent model / main model provided in config. Automatically tracks and adds to workspace.")
    parser.add_argument("-i", dest = 'config', help = 'Config file')
    parser.add_argument("--f", dest = 'force', action = "store_true", help = 'whether to recompute the centerlines')
    args = parser.parse_args()
    
    # manager
    M = Manager(args.config)
    
    # force unregister things
    if args.force:
        # unregister
        c = M.unregister('centerlines', depth=['workspace'])
        o = M.unregister('outlet_file', depth=['workspace'])
        if c:
            Path(c).unlink()
        if o:
            Path(o).unlink()

    
    # retrieve files
    mdl = M['workspace']['mdl']
    vtp = M['workspace']['surface_model']
    inlet = M['metadata']['inlet']
    
    # add centerlines to manager
    cent_file = Path(M['workspace']['root']) / (M['metadata']['model_name'] + '_centerlines.vtp')
    if cent_file.exists():
        print("Centerline files already exist. Run using --f flag to recompute.")
        exit(1)
    M.register(key='centerlines', value=str(cent_file), depth=['workspace'])
    
    # construct centerlines
    c = Centerlines(None)
    ordered_outlets = c.generate_centerlines(mdl, vtp, inlet, str(cent_file))
    
    # write outlets and add to manager
    outlets = str(M.root / "outlet_face_names.dat")
    with open(outlets, "w") as ofile:
        for o in ordered_outlets:
            ofile.write(str(o) + '\n')
    M.register(key='outlet_file', value=outlets, depth=['workspace'])
    
    # update manager
    M.update()
    