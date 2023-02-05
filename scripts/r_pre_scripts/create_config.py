# File: create_config.py
# File Created: Tuesday, 1st November 2022 1:18:06 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 23rd January 2023 7:40:16 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Creates a template config.json in the file with the initial format for the pipeline. Ran before to ensure all files exist in the pipeline.

from sgt.core.manager import Manager
from sgt.utils.parser import Parser
from pathlib import Path

if __name__ == '__main__':
    
    parser = Parser(desc = 'Constructs a template config.json file in directory')
    parser.parser.add_argument('-i', dest = 'dir', help = 'directory to write config file to.')
    
    args = parser.parse_args()
    
    Manager.construct_dev_config(Path(args.dir))
    