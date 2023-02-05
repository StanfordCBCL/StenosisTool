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