# File: setup_workspace.py
# File Created: Monday, 31st October 2022 6:47:53 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Monday, 13th February 2023 2:15:31 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Sets up a results workspace from a SV Model Directory


from sgt.utils.parser import ToolParser
from sgt.core.manager import SVManager

from pathlib import Path

if __name__ == '__main__':
    
    
    
    parser = ToolParser(desc = 'Sets up workspace from SV Model Directory')
    parser.parser.add_argument('-o', dest = 'outdir', help = 'Output directory folder to start workspace.')
    parser.parser.add_argument('-f', dest = 'force', action = 'store_true', default = False, help = 'force re-set up workspace')
    
    args = parser.parse_args()
    
    sv = SVManager(args.config)
    sv.construct_workspace(Path(args.outdir), args.force)
    print('Successfully constructed workspace:', args.outdir)
    
    