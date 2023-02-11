
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
    
    def _update(self):
        """Private method to update yaml file.
        """
        self.write(self.yaml_file)
    
    def __repr__(self) -> str:
        return str(self.yaml)
    
    ###############################
    # Config Manipulation Methods #
    ###############################
    
    def __getitem__(self, key: str):
        """Retrieve a valid registered filepath

        Args:
            key (str): name of registered filepath

        Raises:
            ValueError: Value is not a string.

        Returns:
            str: requested filepath
        """
        # try to find key
        if key not in self.yaml:
            return None
        # key found
        rez = self.yaml[key]
        if not rez:
            # key empty
            return None
        elif type(rez) == str:
            # valid filepath
            return str(self.root / rez)
        else:
            # somehow reached invalid value
            raise ValueError(f"{rez} should be a string.")
    
    def __setitem__(self, key: str, value: Union[Path, str]):
        """registers an item in yaml.

        Args:
            key (str): key to set
            value (Union[Path, str]): file path to insert

        Raises:
            ValueError: value is not a string or Path object
        """

        if type(value) in [Path, str]: 
            # feeding a str -> new filepath
            self.yaml[key] = str(value)
        else:
            raise ValueError("Value must be of Path or str.")
        self._update()
    
    def unregister(self, key: str):
        """unregisters an item from yaml

        Args:
            key (str): key to unregister

        Returns:
            str: filepath of unregistered key
        """
        if key in self.yaml:
            tmp = self.yaml.pop(key)
            self._update()
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
        