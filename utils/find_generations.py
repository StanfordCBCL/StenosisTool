from collections import deque, defaultdict

from .file_io import Solver0D


def find_generations(solver0d: Solver0D):
    ''' Takes a solver0d file and constructs a dictionary mapping of which branches are in which generation '''
    
    inflow_vess, _ = solver0d.identify_inflow_vessel()
    junction_matrix = solver0d.junc_mat
    vess_map = solver0d.vessel_id_map
    num_vess = len(vess_map)
        
    
    # solve for initial cond (gen 0)
    init_list = [vess_map[inflow_vess]]
    temp_vess = inflow_vess
    while len(junction_matrix[temp_vess]) == 1:
        init_list.append(vess_map[temp_vess])
        temp_vess = junction_matrix[temp_vess][0]
        
    generation_map = defaultdict(dict)
    generation_map[0][0] = init_list

    # use a BFS search through junction matrix to find all generations
    cur_gen = 1
    cur_br = 0
    bfs_queue = deque()
    bfs_queue.append(temp_vess)
    bfs_queue.append('end_gen')
    while bfs_queue: 
        cur_vessid = bfs_queue.popleft()
        if cur_vessid == 'end_gen':
            cur_gen += 1
            cur_br = 0
            if bfs_queue:
                bfs_queue.append('end_gen')
            continue
        
        
        
        for vessid in junction_matrix[cur_vessid]:
            segment_list = []
            segment_list.append(vess_map[vessid])
            while len(junction_matrix[vessid]) == 1: # if it is an internal junction, continue until a normal junction  or end is reached.
                vessid = junction_matrix[vessid][0]
                segment_list.append(vess_map[vessid])
            generation_map[cur_gen][cur_br] = segment_list
            bfs_queue.append(vessid)
            cur_br += 1
            
    mapped_num_vess = 0
    for gen, branches in generation_map.items():
        for bid, br in branches.items():
            for vess in br:
                mapped_num_vess+=1
    assert mapped_num_vess == num_vess, 'Number of mapped vessels does not equal number of total vessels'
    print(mapped_num_vess, 'vessels were mapped')
        
    return generation_map

def index_most_stenosed(vessel_list: list)-> int:
    ''' in a list of vessel dicts, finds index with greatest stenosis factor'''
    max_index = 0
    cur_index = 1
    cur_max = vessel_list[0]['zero_d_element_values']['stenosis_coefficient']
    for vess in vessel_list[1:]:
        if vess['zero_d_element_values']['stenosis_coefficient'] > cur_max:
            cur_max = vess['zero_d_element_values']['stenosis_coefficient']
            max_index = cur_index
        cur_index+=1
    return max_index

    

    
    