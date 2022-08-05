import os
import argparse
import sys

def m2d(val):
    '''convert mmHg to dynes/cm^2
    '''
    return val * 1333.22

def d2m(val):
    '''convert dynes/cm^2 to mmHg
    '''
    return val / 1333.22

def blockPrint():
    ''' block printing 
    '''
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    ''' enable printing
    '''
    sys.stdout = sys.__stdout__
    
def create_parser(desc = '') -> argparse.ArgumentParser:
    ''' creates a parse for scripts '''
    parser = argparse.ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(help='Mode for script', dest = 'mode')
    tool = subparsers.add_parser(name='tool')
    dev = subparsers.add_parser(name='dev')
    dev.add_argument('-root', default='.', help='Root of entire project')
    dev.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    return parser, dev, tool

def create_dev_parser(desc = '') -> argparse.ArgumentParser:
    ''' creates a parse for scripts '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-root', default='.', help='Root of entire project')
    parser.add_argument('-models', dest = 'models', nargs = '*', default = [], help = 'Specific models to run')
    return parser

def create_tool_parser(desc = '') -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-solver_dirs', default=[], nargs = '*', help='Solver Directories to begin searching from')
    return parser

def get_solver_name(solver_dir):
    for name in os.listdir(solver_dir):
        if os.path.splitext(name)[-1] == '.in':
            return name
        
def get_solver_path(solver_dir):
    solver_name = get_solver_name(solver_dir)
    if solver_name is None:
        return
    return os.path.join(solver_dir, solver_name)


def get_basename(fp):
    return os.path.splitext(os.path.basename(fp))[0]