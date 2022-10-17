
import argparse

class Parser():
    ''' base class for our parsers
    '''
    def __init__(self, desc):
        self.parser = argparse.ArgumentParser(description=desc)
    
    def parse_args(self):
        return self.parser.parse_args()

class ToolParser(Parser):
    ''' parser for tool-based scripts 
    '''
    def __init__(self, desc):
        super().__init__(desc)
        self.parser.add_argument('-solver_dirs', default=[], nargs = '*', help='Solver Directories to begin searching from')

class DevParser(Parser):
    ''' parser for dev scripts
    '''
    def __init__(self, desc):
        super().__init__(desc)
        self.parser.add_argument('-root', default='.', help='Root of entire project')
        self.parser.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')