import os
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
    