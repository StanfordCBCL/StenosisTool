# File: generate_sv_config.py
# File Created: Monday, 13th February 2023 2:27:57 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Thursday, 13th July 2023 1:37:31 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Generates a template config for generating a workspace using setup_workspace.py.

import sys
from pathlib import Path
import shutil
import argparse

if __name__ == '__main__':
  
    parser = argparse.ArgumentParser(description="Generates a template config for generating a workspace using setup_workspace.py.")
    parser.add_argument("-o", dest="outfile", required=True, 
                        help="The path to where the config file should go. Specify a valid filename to rename the file, otherwise a default 'config.yaml' is used.")
    args = parser.parse_args()
    
    # writes config file
    shutil.copy(str(Path(__file__).parent / "config_template" / "config.yaml"), args.outfile)