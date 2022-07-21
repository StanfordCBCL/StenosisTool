import os
import configparser
from collections import defaultdict, deque
import xml.etree.ElementTree as ET
import json
import numpy as np

class Solver0D():
    ''' Solver File modification '''
    
    BC = 'boundary_conditions'
    VESS = 'vessels'
    SIM = 'simulation_parameters'
    JUNC = 'junctions'
    DESC = 'description'
    
    class Node():
        
        def __init__(self, vess_id, vessel_info, generation = 0):
            self.vess_id = vess_id
            self.vessel_info = vessel_info
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
    
    def get_vessel_tree(self):
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
        ''' tree to list of vessel ids'''

        vessel_list = []
        for node in self.tree_bfs_iterator(tree):
            vessel_list.append(node)
                
        return vessel_list      
            
            
        
        
        
        
        
    
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
    ''' 0D solver results file '''
    def __init__(self, results_file):
        self.results = np.load(results_file, allow_pickle= True).item()
        self.pressures = self.results['pressure']
        self.time = self.results['time']
        self.flow = self.results['flow']



class NotFound(Exception):
    pass

def check_exists(file, err = None, mkdir = False):
    ''' check if file/path exists'''
    if err == None:
        err = file + ' does not exist'
    if not os.path.exists(file):
        if mkdir:
            print(err + ': making directory')
            os.mkdir(file)
        else:
            raise NotFound(err)
    return file

def check_exists_bool(file, err = None, ignore = False):
    ''' check if file/path exists'''
    if err == None:
        err = file + ' does not exist'
    if not os.path.exists(file):
        if ignore:
            print(err)
            return False
        else:
            raise NotFound(err)
    return True


def parse_config(config_file, default_section = configparser.DEFAULTSECT):
    ''' turn a config file to a dict (all strs) '''
    config = configparser.ConfigParser(allow_no_value=True,
                                       default_section=default_section,
                                       interpolation=configparser.ExtendedInterpolation())
    config.read(config_file)
    
    cfg_dict = defaultdict(dict)
    for section in config.sections():
        for option in config.options(section):
            cfg_dict[section][option] = config.get(section, option)
    return cfg_dict

def parse_mdl(mdl_file, reverse = False):
    ''' parses mdl file for faceid corresponding to caps'''
    
    # since .mdl has a line like "<format version="1.0" />" which fails for the standard XMLparser, rather than creating a custom parser, just remove that line after reading the file in and parse as a list of strings
    with open(mdl_file, 'r') as mdl:
        lines = mdl.readlines()
        if 'format' in lines[1]:
            lines[1] = ''
    
    root = ET.fromstringlist(lines)
    faces = root[0][0].find('faces')

    # save to a dict
    face_mappings = {}
    for face in faces:
        if face.attrib['type'] == 'cap':
            if reverse:
                face_mappings[int(face.attrib['id'])] = face.attrib['name']
            else:
                face_mappings[face.attrib['name']] = int(face.attrib['id'])
    return face_mappings

def write_json(fp, data):
    with open(fp, 'w') as sfile:
            json.dump(data, sfile, indent = 4, sort_keys=True)