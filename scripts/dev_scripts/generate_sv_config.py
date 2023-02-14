# File: generate_sv_config.py
# File Created: Monday, 13th February 2023 2:27:57 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 13th February 2023 2:30:54 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Generates initial config from the Simvascular Folder

import sys
from pathlib import Path

yaml = """# SV Workspace
workspace:
  surface_model: "Models/0080_0001.vtp"
  centerlines: "centerline_files/0080_0001_model_centerlines.vtp"
  mdl: "Models/0080_0001.mdl"
  flow_file: "flow-files/inflow_1d.flow"
  rcrt_file: 
  capinfo: "Models/CapInfo"

# Metadata about model
metadata:
  model_name: "0080_0001"
  diseased: false
  inlet: "inflow"

# Options in pipeline
options:
  tune: true

# Tuning Parameters
tune_params:
  maxPAP:
    - 18 
    - 25
  minPAP: 
    - 8
    - 12
  mPAP: 
    - 12
    - 16
  rpa_split: 0.55
  PCWP: 7"""

if __name__ == '__main__':
    
    outfile = Path(sys.argv[1])
    # writes config file at sys.argv[1]
    with outfile.open("w") as ofp:
        ofp.write(yaml)