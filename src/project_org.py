
from collections import defaultdict
from os import path
import os

from .misc import *

class ProjectPath():
    ''' Entire Project Manager'''
    
    class DataPath():
        ''' Path manager for data folder'''
        
        HEALTHY = 'healthy'
        NCI_STEN = 'nci_stenosis'
        CI_STEN = 'ci_stenosis'
        
        
        
        class ModelPath():
            ''' Model specific information should be configured under config_files in dev_config.ini'''
            
            DEV_CONF = 'dev_config.ini'
            CONFIG = 'config_files'
            CENTERLINES = 'centerline_files'
            SOLVER = 'solver_files'
            
            MOD_CENT = 'model_centerlines.vtp'
            TUNE_CENT = 'tune_centerlines.vtp'
            
            MODEL_SOLVER = 'model.in'
            TUNE_SOLVER = 'tuner.in'
            
            ARTIFICIAL_STEN = 'artificial_stenosis'
            
            STENOSIS_FILE = 'stenosis_file.txt'
            STENOSIS_SANITY_DIR = 'stenosis_sanity_dir'
            
            
            def __init__(self, root, type):
                self.model_root = root
                
                self.config_files = check_exists(path.join(root, self.CONFIG), 
                                                 err = 'Model {} does not have a {} directory'.format(path.basename(root), self.CONFIG))
                
                # search and read for dev_config file in self.config_files
                self.dev_config = check_exists(path.join(self.config_files, self.DEV_CONF), 
                                               err = 'Model {} does not have a {} file. Please create one.'.format(os.path.basename(root), self.DEV_CONF))
                
                # read dev config
                self.params = read_config(self.dev_config)
                
                for key, vals in self.params['files'].items():
                    if vals != '':
                        self.params['files'][key] = path.join(root, vals)
                                    
                self.centerline_files = check_exists(path.join(root, self.CENTERLINES), mkdir = True)
                self.solver_files = check_exists(path.join(root, self.SOLVER), mkdir = True)
                
                self.model_solver = path.join(self.solver_files, self.params['metadata']['name'] + '_' + self.MODEL_SOLVER)
                self.tune_solver = path.join(self.solver_files, self.params['metadata']['name'] + '_' + self.TUNE_SOLVER)
                self.model_centerlines = path.join(self.centerline_files, self.params['metadata']['name'] + '_' + self.MOD_CENT)
                self.tune_centerlines = path.join(self.centerline_files, self.params['metadata']['name'] + '_' + self.TUNE_CENT)
                
                self.type = type
                if type == 'healthy':
                    self.artificial_stenosis_dir = check_exists(path.join(root, self.ARTIFICIAL_STEN), mkdir = True)
                    self.stenosis_vessels_file = None
                elif type == 'nci_stenosis':
                    self.artificial_stenosis_dir = None
                    self.stenosis_vessels_file = path.join(self.solver_files, self.STENOSIS_FILE)
                    self.stenosis_check_dir = check_exists(path.join(root, self.STENOSIS_SANITY_DIR), mkdir = True)
                elif type == 'ci_stenosis':
                    self.artificial_stenosis_dir = None
                    self.stenosis_vessels_file = path.join(self.solver_files, self.STENOSIS_FILE)
                
            def __repr__(self):
                s = ''
                for attr, val in self.__dict__.items():
                    s += '{} : \n\t{}\n'.format(attr, val)
                return s
        

        
        def __init__(self, root):
            
            self.data_root = root
            self.healthy = check_exists(path.join(root, self.HEALTHY), mkdir = True)
            self.nci_stenosis = check_exists(path.join(root, self.NCI_STEN), mkdir = True)
            self.ci_stenosis = None # ! check_exists(path.join(root, self.CI_STEN), mkdir = True)
            self.data_dirs = [self.healthy, self.nci_stenosis, self.ci_stenosis]
            
            self.models = defaultdict(dict)
            for idx, (name, dir) in enumerate(zip([self.HEALTHY, self.NCI_STEN, self.CI_STEN], self.data_dirs)):
                if dir == None:
                    continue
                for mod in os.listdir(dir):
                    model_path = path.join(dir, mod)
                    if os.path.isdir(model_path):
                        self.models[name][mod] = self.ModelPath(model_path, type = name)
        
        def __repr__(self):
            s = ''
            for attr, val in self.__dict__.items():
                s += '{} : \n\t{}\n'.format(attr, val)
            return s[1:]

    

        
    class TestPath():
        pass
    
    def __init__(self, root = '.') -> None:
        
        self.data_dir = check_exists(path.join(root, 'data'), mkdir = True)
        self.data_path = self.DataPath(self.data_dir)
        
        self.test_dir = check_exists(path.join(root, 'tests'), mkdir = True)
        
        self.logs_dir = check_exists(path.join(root, 'logs'), mkdir = True)
        self.female_bsa_chart = parse_bsa_chart(path.join(self.logs_dir, 'female_bsa_chart.dat'))
        self.male_bsa_chart = parse_bsa_chart(path.join(self.logs_dir, 'male_bsa_chart.dat'))
        self.female_height_chart = parse_height_chart(path.join(self.logs_dir, 'female_heights.dat'))
        self.male_height_chart = parse_height_chart(path.join(self.logs_dir, 'male_heights.dat'))
        
        
        
        
    def run_all(self, func, *args):
        ''' wrapper for functions utilizing all models. Function must reserve first position for a ModelPath object'''
        output = {}
        for cond, modeldict in self.data_path.models.items():
            for mod_name, mod_path in modeldict.items():
                print('Running function for {}...'.format(mod_name))
                output[mod_name] = func(mod_path, *args)
        return output
    
    def find_model(self, model_name) -> DataPath.ModelPath:
        if model_name in self.data_path.models[self.data_path.HEALTHY]:
            return self.data_path.models[self.data_path.HEALTHY][model_name]
        elif model_name in self.data_path.models[self.data_path.CI_STEN]:
            return self.data_path.models[self.data_path.CI_STEN][model_name]
        elif model_name in self.data_path.models[self.data_path.NCI_STEN]:
            return self.data_path.models[self.data_path.NCI_STEN][model_name]
        else:
            print(model_name + ' is an invalid model')
            return None
        
            
    
    def __repr__(self):
        s = ''
        for attr, val in self.__dict__.items():
            s += '{} : \n\t{}\n'.format(attr, val)
        return s
        
            
        

        
if __name__ == '__main__':
    import sys
    test = ProjectPath(sys.argv[1])
    
    print(test)