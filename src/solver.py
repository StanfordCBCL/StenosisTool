import numpy as np
import json
from collections import defaultdict, deque
import pandas as pd

from .file_io import write_json


class Solver0D():
    ''' Solver File modification '''
    
    BC = 'boundary_conditions'
    VESS = 'vessels'
    SIM = 'simulation_parameters'
    JUNC = 'junctions'
    DESC = 'description'
    
    class Node():
        ''' node for vessel tree '''
        def __init__(self, vess_id, vessel_info, generation = 0):
            self.vess_id = vess_id
            self.vessel_info = vessel_info
            self.parent = None
            self.children = None
            self.generation = generation

    class BranchNode():
        ''' node for branch tree'''
        def __init__(self, vess_id, vessel_info, generation):
            self.branch_id = int(vessel_info['vessel_name'].split('_')[0][6:])
            self.vess_id = [vess_id]
            self.vessel_info = [vessel_info]
            self.parent = None
            self.children = None
            self.generation = generation
    
    
    def __init__(self):
        self.solver_file = None
        self.solver_data = None
        self._vessel = None
        self._bc = None
        self._description = None
        self._simulation_params = None
        self._junctions = None
        self.junc_mat = None
        self.vessel_id_map = None
        self.vessel_name_map = None

        
    def set_up_new_solver(self):
        ''' sets up a new solver'''
        self.solver_data = {self.BC: [], 
                            self.VESS: [],
                            self.SIM: {},
                            self.JUNC: [],
                            self.DESC: {}}
        self.update_solver_data()
        
    
    def read_solver_file(self, solver_file):
        ''' reads the solver file into a dict'''
        self.solver_file = solver_file
        with open(solver_file, 'r') as sfile:
            self.solver_data = json.load(sfile)
        self.update_solver_data()
    
    def update_solver_data(self):
        ''' update attributes if solver data changes'''
        self._vessel = self.solver_data[self.VESS]
        self._bc = self.solver_data[self.BC]
        if self.DESC in self.solver_data:
            self._description = self.solver_data[self.DESC]
        else:
            self._description = {}
        self._simulation_params = self.solver_data[self.SIM]
        self._junctions = self.solver_data[self.JUNC]
        if self._bc:
            for bc in self._bc:
                if bc['bc_type'] == 'FLOW':
                    break
            self.tc = bc['bc_values']['t'][-1] - bc['bc_values']['t'][0]
        else:
            self.tc = 0
        self._construct_vessel_map()
        self._junction_matrix()

    def write_solver_file(self, solver_file):
        ''' writes a dict into the solver file'''
        write_json(solver_file, self.solver_data)
        
    
    # Contruct a vessel Map to avoid linear search
    def _construct_vessel_map(self):
        self.vessel_id_map = {}
        self.vessel_name_map = {}
        for vess in self.vessel:
            self.vessel_id_map[vess['vessel_id']] = vess
            self.vessel_name_map[vess['vessel_name']] = vess

        
    # retrieve a modifiable dict of a vessel
    def get_vessel(self, id = None, name = None):
        ''' retrieves a vessel dict that can be modified directly'''
        if id is None and name is None:
            print('No identification was specified')
            return None
        if id is not None:
            return self.vessel_id_map[id]
        if name is not None:
            return self.vessel_name_map[name]
        
        
    # Obtain a node matrix of junction to vessel id
    def _junction_matrix(self):
        self.junc_mat = defaultdict(list)
        for junc in self.junctions:
            self.junc_mat[junc['inlet_vessels'][0]] = junc['outlet_vessels']
        
    def identify_inflow_vessel(self):
        for vess in self.vessel:
            if 'boundary_conditions' in vess:
                if 'inlet' in vess['boundary_conditions']:
                    return vess['vessel_id'], vess['vessel_name']
        print('Error: an inlet bc does not exist.')
        return 
    
    def get_vessel_tree(self) -> Node:
        ''' get the vessel tree at that snapshot in the class '''
        
        head_node = self.Node(0, self.get_vessel(0), generation=0)
        head_node.children = [self.junc_mat[0]]
        
        
        bfs_queue = deque()
        bfs_queue.append(head_node)
        
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            if len(self.junc_mat[cur_node.vess_id]) != 1: # if not an internal junction
                next_gen = cur_node.generation + 1
            else:
                next_gen = cur_node.generation # stays in current generation if it is an internal junction
                
            cur_node.children = [self.Node(child_id, self.get_vessel(child_id), generation = next_gen) for child_id in self.junc_mat[cur_node.vess_id]]
            for child_node in cur_node.children:
                child_node.parent = cur_node.vess_id
                bfs_queue.append(child_node)
        
        return head_node
    
    def get_branch_tree(self) -> BranchNode:
        ''' get the branch tree at that snapshot in the class '''
        
        head_node = self.BranchNode(vess_id=0, 
                                    vessel_info=self.get_vessel(0), 
                                    generation=0)
        
        
        bfs_queue = deque()
        bfs_queue.append(head_node)
    
        
        while bfs_queue:
            cur_node = bfs_queue.popleft()
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
        
        return head_node
    
    def tree_bfs_iterator(self, tree):
        ''' iterates using bfs'''
        bfs_queue = deque()
        bfs_queue.append(tree)
        while bfs_queue:
            cur_node = bfs_queue.popleft()
            for child_node in cur_node.children:
                bfs_queue.append(child_node)
            yield cur_node
    
    def tree_to_list(self, tree):
        ''' tree to list of vessel nodes'''

        vessel_list = []
        for node in self.tree_bfs_iterator(tree):
            vessel_list.append(node)
                
        return vessel_list      
            
            
    def to_cpp(self):
        for junc in self.junctions:
            if junc['junction_type'] == 'BloodVesselJunction':
                junc['junction_type'] = 'NORMAL_JUNCTION'
        
        
    
    def num_vessel_segments(self):
        return len(self.vessel)
    
    def num_junctions(self):
        return len(self.junctions)
    
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
        
        
class SolverResults():
    ''' 0D c solver results file '''
    def __init__(self, df: pd.DataFrame, last_cycle = False, tc = 0):
        
        if last_cycle:
            # only save last cycle
            df = self.only_last_cycle(df, tc)
        self.result_df = df
        
    @staticmethod
    def only_last_cycle(df, tc) -> pd.DataFrame:
        ''' convert to only last cycle '''
        df = df.copy()
        df = df[df['time'] >= df['time'].max() - tc].copy()
        df['time'] -= df['time'].min()
        return df
       
    @classmethod
    def load_from_csv(cls, fp):
        df = pd.read_csv(fp)
        return cls(df, last_cycle=False, tc = 0)
    
    def vessel_df(self, vessel_name):
        ''' retrieves a df isolated by name'''
        return self.result_df[self.result_df['name'] == vessel_name]
    
    def convert_to_branch_results(self, solver: Solver0D, out_file):
        ''' Convert the c results into branch results (not sorted)'''
        branch_results = {}
        
        # write time
        branch_results['time'] = np.array(list(self.vessel_df('V0')['time']))
        
        #! Not entering them sorted (shouldnt matter)
        branch_results['distance'] = {}
        branch_results['flow'] = {}
        branch_results['pressure'] = {}
        
        branch_tree = solver.get_branch_tree()

        for node in solver.tree_bfs_iterator(branch_tree):
            branch_id = node.branch_id
            # get initial values for first vessel
            vess_Q = []
            vess_P = []
            vess_D = [0.0]
            vess_df = self.vessel_df('V' + str(node.vess_id[0]))
            vess_Q.append(list(vess_df['flow_in']))
            vess_Q.append(list(vess_df['flow_out']))
            vess_P.append(list(vess_df['pressure_in']))
            vess_P.append(list(vess_df['pressure_out']))
            vess_D.append(node.vessel_info[0]['vessel_length'])
            
            # if longer than just 1 vessel segment per branch
            if len(node.vess_id) > 1:
                cur_idx = 1
                for vess_id in node.vess_id[1:]:
                    vess_df = self.vessel_df('V' + str(vess_id))
                    vess_Q.append(list(vess_df['flow_out']))
                    vess_P.append(list(vess_df['pressure_out']))
                    vess_D.append(vess_D[-1] + node.vessel_info[cur_idx]['vessel_length'])
                    cur_idx += 1
        
            branch_results['distance'][branch_id] = np.array(vess_D)
            branch_results['flow'][branch_id] = np.array(vess_Q)
            branch_results['pressure'][branch_id] = np.array(vess_P)
        
    
        np.save(out_file, branch_results, allow_pickle=True)
        
    def save_csv(self, out_file):
        self.result_df.to_csv(out_file, sep = ',', header = True, index = False)
        

