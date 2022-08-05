
import json
from collections import defaultdict, deque
import numpy as np

from .file_io import write_json
from .flow import Inflow0D
from .bc import BoundaryConditions
        
        
class Solver0D():
    ''' Solver File modification '''
    
    # config file params
    BC = 'boundary_conditions'
    VESS = 'vessels'
    SIM = 'simulation_parameters'
    JUNC = 'junctions'
    DESC = 'description'
    BC_MAP = 'bc_map'
    
    class Node():
        ''' base node for trees
        '''
        def __init__(self):
            self.vess_id = []
            self.vessel_info = []
            self.parent = None
            self.children = []
            self.generation = -1
            self.side = None
        
        def last_vessel(self):
            ''' get the last vessel '''
            return self.vess_id[-1], self.vessel_info[-1]
        
        def get_type(self):
            ''' if Node is an outlet '''
            lvid, lvinfo = self.last_vessel()
            if 'boundary_conditions' in lvinfo:
                if 'inlet' in lvinfo['boundary_conditions']:
                    return 'inlet'
                if 'outlet' in lvinfo['boundary_conditions']:
                    return 'outlet'
            return 'internal'
        
    class VesselNode(Node):
        ''' node for vessel tree
        '''
        def __init__(self, vess_id, vessel_info, generation):
            super().__init__()
            self.vess_id.append(vess_id)
            self.vessel_info.append(vessel_info)
            self.generation = generation

    class BranchNode(Node):
        ''' node for branch tree
        '''
        def __init__(self, vess_id, vessel_info, generation):
            super().__init__()
            self.vess_id.append(vess_id)
            self.vessel_info.append(vessel_info)
            self.branch_id = self._branch_id()
            self.generation = generation
            
        def _branch_id(self):
            return int(self.vessel_info[-1]['vessel_name'].split('_')[0][6:])
        
        def get_branch_len(self):
            length = 0
            for vessel in self.vessel_info:
                length += vessel['vessel_length']
            return length
    
                
    def __init__(self):
        self.solver_file = None
        self.solver_data = None
        self._vessel = None
        self._bc = None
        self._bc_map = None
        self._description = None
        self._simulation_params = None
        self._junctions = None
        self.inflow = None
        self.junc_mat = None
        self.vessel_id_map = None
        self.vessel_name_map = None

    # Solver IO
    def setup_empty_solver(self):
        ''' sets up a new solver
        '''
        self.solver_data = {self.BC: [],
                            self.BC_MAP: {}, 
                            self.VESS: [],
                            self.SIM: {},
                            self.JUNC: [],
                            self.DESC: {}}
        self._update_solver_data()
        
    def write_solver_file(self, solver_file):
        ''' writes a dict into the solver file'''
        write_json(solver_file, self.solver_data)
    
    def read_solver_file(self, solver_file):
        ''' reads the solver file into a dict
        '''
        self.solver_file = solver_file
        with open(solver_file, 'r') as sfile:
            self.solver_data = json.load(sfile)
        self._update_solver_data()
    
    # Private Operations
    def _update_solver_data(self):
        ''' update attributes if solver data changes
        '''
        # map to attributes
        self._vessel = self.solver_data[self.VESS]
        self._bc = self.solver_data[self.BC]
        if self.DESC in self.solver_data:
            self._description = self.solver_data[self.DESC]
        else:
            self._description = {}
        self._simulation_params = self.solver_data[self.SIM]
        self._junctions = self.solver_data[self.JUNC]
        if self.BC_MAP in self.solver_data:
            self._bc_map = self.solver_data[self.BC_MAP]
        else:
            self._bc_map = {}
        
        # compute inflow if it exists. Otherwise use a dummy.
        for bc in self._bc:
            if bc['bc_type'] == 'FLOW':
                self.inflow = Inflow0D(inflow_arr = np.array(list(zip(bc['bc_values']['t'], bc['bc_values']['Q']))))
                break
        if self.inflow is None:
            self.inflow = Inflow0D(inflow_arr = np.array([(0,0),(1,0)]), smooth = False)
            
        # construct vessel maps and junction matrix
        self._construct_vessel_map()
        self._junction_matrix()

    def _construct_vessel_map(self):
        ''' Vessel Map to avoid Linear searches
        '''
        self.vessel_id_map = {}
        self.vessel_name_map = {}
        for vess in self.vessel:
            self.vessel_id_map[vess['vessel_id']] = vess
            self.vessel_name_map[vess['vessel_name']] = vess
            
    def _junction_matrix(self):
        ''' Junction Maxtrix to map inlet and outlet of junctions'''
        self.junc_mat = defaultdict(list)
        for junc in self.junctions:
            self.junc_mat[junc['inlet_vessels'][0]] = junc['outlet_vessels']

    def _det_lpa_rpa(self, head_node):
        
        if self.bc_map is None:
            return
        
        head_node.side = 'mpa'
        side1 = head_node.children[0]
        side2 = head_node.children[1]
        
        cur_branch = side1
        while cur_branch.children:
            cur_branch = cur_branch.children[0]
    
        rcr_name = self.bc_map[cur_branch.vessel_info[-1]['boundary_conditions']['outlet']]
        if 'lpa' in rcr_name.lower():
            for node in self.tree_bfs_iterator(side1):
                node.side = 'lpa'
            for node in self.tree_bfs_iterator(side2):
                node.side = 'rpa'
        elif 'rpa' in rcr_name.lower():
            for node in self.tree_bfs_iterator(side1):
                node.side = 'rpa'
            for node in self.tree_bfs_iterator(side2):
                node.side = 'lpa'
        else:
            raise ValueError(rcr_name + ' is an invalid rcr name')

    # Public Functions
    
    def get_vessel(self, id = None, name = None):
        ''' retrieves a vessel dict that can be modified directly'''
        if id is None and name is None:
            print('No identification was specified')
            return None
        if id is not None:
            return self.vessel_id_map[id]
        if name is not None:
            return self.vessel_name_map[name]
        
    def identify_inflow_vessel(self):
        ''' identify the inlet vessel'''
        for vess in self.vessel:
            if 'boundary_conditions' in vess:
                if 'inlet' in vess['boundary_conditions']:
                    return vess['vessel_id'], vess['vessel_name']
        print('Error: an inlet bc does not exist.')
        return None
    
    def get_vessel_tree(self) -> VesselNode:
        ''' get a vessel tree at that snapshot in the class
        '''
        
        head_node = self.VesselNode(0, self.get_vessel(0), generation=0)
        
        bfs_queue = deque()
        bfs_queue.append(head_node)
        
        # bfs traverse
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            if len(self.junc_mat[cur_node.vess_id[0]]) != 1: # if not an internal junction
                next_gen = cur_node.generation + 1
            else:
                next_gen = cur_node.generation # stays in current generation if it is an internal junction
                
            cur_node.children = [self.VesselNode(child_id, self.get_vessel(child_id), generation = next_gen) for child_id in self.junc_mat[cur_node.vess_id[0]]]
            for child_node in cur_node.children:
                child_node.parent = cur_node.vess_id[0]
                bfs_queue.append(child_node)
        
        # determine which side each node is on
        self._det_lpa_rpa(head_node)
        
        return head_node
    
    def get_branch_tree(self) -> BranchNode:
        ''' get the branch tree at that snapshot in the class '''
        
        head_node = self.BranchNode(vess_id=0, 
                                    vessel_info=self.get_vessel(0), 
                                    generation=0)
        
        bfs_queue = deque()
        bfs_queue.append(head_node)
    
        # bfs traverse
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            # scan for internal junctions
            while True:
                # junc_mat
                child_vess = self.junc_mat[cur_node.vess_id[-1]]
                # internal junction so add to same node
                if len(child_vess) == 1:
                    cur_node.vess_id.append(child_vess[0])
                    cur_node.vessel_info.append(self.get_vessel(child_vess[0]))
                # otherwise, exit while loop
                else:
                    break
                
            # determine children
            cur_node.children = [self.BranchNode(vess_id=child_id, vessel_info=self.get_vessel(child_id), generation=cur_node.generation + 1) for child_id in child_vess]
            for child_node in cur_node.children:
                child_node.parent = cur_node
                bfs_queue.append(child_node)
        
        self._det_lpa_rpa(head_node)
        
        return head_node
    
    def tree_bfs_iterator(self, tree) -> Node:
        ''' Tree BFS Iterator 
        '''
        bfs_queue = deque()
        bfs_queue.append(tree)
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            for child_node in cur_node.children:
                bfs_queue.append(child_node)
            yield cur_node
    
    def tree_to_list(self, tree):
        ''' Tree to list of vessel nodes
        '''

        vessel_list = []
        for node in self.tree_bfs_iterator(tree):
            vessel_list.append(node)
                
        return vessel_list      
            
    def to_cpp(self):
        ''' convert to cpp solver structure
        '''
        for junc in self.junctions:
            if junc['junction_type'] == 'BloodVesselJunction':
                junc['junction_type'] = 'NORMAL_JUNCTION'
    
    def write_bc_map(self, boundary_conditions: BoundaryConditions):
        ''' must be in same order as bc that was written to rcrt.dat file
        '''
        bc_map = {}
        counter = 0
        for bc in boundary_conditions.bc_list:
            bc_map[bc['type'] + '_' + str(counter)] = bc['faceID']
            counter+=1
        self.bc_map = bc_map
    
    def num_vessel_segments(self):
        return len(self.vessel)
    
    def num_junctions(self):
        return len(self.junctions)
    
    # Properties
    
    @property
    def vessel(self):
        return self._vessel

    @vessel.setter
    def vessel(self, val):
        self._vessel = val
        self.solver_data[self.VESS] = val
        
    @property
    def bc(self):
        return self._bc

    @bc.setter
    def bc(self, val):
        self._bc = val
        self.solver_data[self.BC] = val
    
    @property
    def simulation_params(self):
        return self._simulation_params

    @simulation_params.setter
    def simulation_params(self, val):
        self._simulation_params = val
        self.solver_data[self.SIM] = val

    @property
    def junctions(self):
        return self._junctions

    @junctions.setter
    def junctions(self, val):
        self._junctions = val
        self.solver_data[self.JUNC] = val
        
    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, val):
        self._description = val
        self.solver_data[self.DESC] = val
    
    @property
    def bc_map(self):
        return self._bc_map

    @bc_map.setter
    def bc_map(self, val):
        self._bc_map = val
        self.solver_data[self.BC_MAP] = val
 
 
 