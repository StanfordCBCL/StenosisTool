# File: generate_AS_config.py
# File Created: Wednesday, 2nd November 2022 2:24:21 pm
# Author: John Lee (jlee88@nd.edu)
# Last Modified: Wednesday, 2nd November 2022 2:47:44 pm
# Modified By: John Lee (jlee88@nd.edu>)
# 
# Description: Generates a config for aritifical stenosis generation, allowing you to select which vessels to stenose.

from sgt.utils.parser import Parser
from sgt.utils.io import write_json

config = {
    'description': ["Inside vessels, please specify which vessel_id and the degree of stenosis as follows:", "[vessel_id, occ_lower, occ_higher] (for an occlusion range)",
    "[vessel_i,d occ] (for a set occlusion)"],
    'random': False,
    'vessels': [
        [0, .75, .9],
        [1, .8]
    ]
}

config_random = {
    'description': ["Inside generations, please specify how much of the length, and to what occlusion range you wish to occlude:", "[generation, percent_length, occ_lower, occ_higher]"],
    'random' : True,
    'generations': [
        [1, .5, .75, .9],
    ]
}



if __name__ == '__main__':
    
    parser = Parser(desc = "Creates a config for artificial stenosis at a specified location ")
    
    parser.parser.add_argument('-o', help = 'output destination path of file')
    parser.parser.add_argument('-r', action = 'store_true', default = False, help = 'Whether to use the random config version')
    
    args = parser.parse_args()
    if args.r:
        write_json(args.o, config_random, sort_keys=False)
    else:
        
        write_json(args.o, config, sort_keys=False)