
import yaml
from pathlib import Path
from typing import Union

class Manager():
    """ A generic manager for file directories that parses a single yaml file. Can register files.
    """

    def __init__(self, yaml_file: str):
        
        # Open and load yaml_file
        self.yaml_file = yaml_file 
        with open(yaml_file, 'r') as yfile:
            self.yaml = yaml.full_load(yfile)
        
        self.root = Path(yaml_file).parent
    
    def write(self, out_file: str):
        """writes yaml file to outfile

        Args:
            out_file (str): output file path + name
        """
        with open(out_file, 'w') as yfile:
            yaml.safe_dump(self.yaml, yfile)
    
    def update(self):
        """Method to update yaml file.
        """
        self.write(self.yaml_file)
    
    def __repr__(self) -> str:
        return str(self.yaml)
    
    ###############################
    # Config Manipulation Methods #
    ###############################
    
    def __getitem__(self, key: str):
        return self.yaml[key]
        
    def register(self, key: str, value, depth: list = []):
        """unregisters an item from yaml

        Args:
            key (str): key to register
            depth (list): depth at which to register

        Returns:
            str: filepath of unregistered key
        """
        # walk down
        final = self.yaml
        for d in depth:
            if d not in final:
                final[d] = {}
            final = final[d]
        
        final[key] = value
    
    def unregister(self, key: str, depth: list = []):
        """unregisters an item from yaml

        Args:
            key (str): key to unregister

        Returns:
            str: filepath of unregistered key
        """
        # walk down
        final = self.yaml
        for d in depth:
            final = final[d]
        #
        if key in final:
            tmp = final.pop(key)
            return tmp
    
    def register_many(self, items: dict):
        """Registers all items in a dictionary

        Args:
            items (dict): a dictionary of key, filepaths
        """
        for key, value in items.items():
            self[key] = value
        
    
    def get_root(self):
        return str(self.root)
    
    ##############
    # IO methods #
    ##############
        