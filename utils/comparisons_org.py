from tabnanny import check
from .project_org import ProjectPath
from .misc import *
from os import path


class ComparisonPath():
    
    OLD_0D = 'old_0d'
    NEW_0D = 'new_0d'
    CONTROL_3D = 'control_3d'
    
    def __init__(self, root = '.') -> None:
        
        self.org = ProjectPath(root)
        self.healthy_pulmonary = self.org.data_path.models['healthy']['0080_0001']
        self.sten_pulmonary_0 = self.org.data_path.models['nci_stenosis']['0118_1000']
        self.sten_pulmonary_1 = self.org.data_path.models['nci_stenosis']['SU0238-AS3']
        
        self.models = [self.healthy_pulmonary, self.sten_pulmonary_0, self.sten_pulmonary_1]
        self.comp_dir = check_exists(path.join(root, 'test_comparisons'), mkdir = True)
        
        self.dir_names = {}
        # set up dirs
        for mod_dir in self.models:
            mod_name = mod_dir.params['metadata']['name']
            self.dir_names[mod_name] = check_exists(path.join(self.comp_dir, mod_name), mkdir = True)
            for i in [self.OLD_0D, self.NEW_0D, self.CONTROL_3D]:
                check_exists(path.join(self.dir_names[mod_name], i), mkdir = True)
            
            
        