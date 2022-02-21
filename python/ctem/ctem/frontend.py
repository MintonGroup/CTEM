import numpy
import os
import subprocess
import shutil
from ctem import io

class Simulation:
    """
    This is a class that defines the basic CTEM simulation object
    """
    def __init__(self, param_file="ctem.in"):
        currentdir = os.getcwd()
        self.parameters = {
            'restart': None,
            'runtype': None,
            'popupconsole': None,
            'saveshaded': None,
            'saverego': None,
            'savepres': None,
            'savetruelist': None,
            'seedn': 1,
            'totalimpacts': 0,
            'ncount': 0,
            'curyear': 0.0,
            'fracdone': 1.0,
            'masstot': 0.0,
            'interval': 0.0,
            'numintervals': 0,
            'pix': -1.0,
            'gridsize': -1,
            'seed': 0,
            'maxcrat': 1.0,
            'shadedminhdefault': 1,
            'shadedmaxhdefault': 1,
            'shadedminh': 0.0,
            'shadedmaxh': 0.0,
            'workingdir': currentdir,
            'ctemfile': param_file,
            'datfile': 'ctem.dat',
            'impfile': None,
            'sfdcompare': None,
            'sfdfile': None
        }
        self.parameters = io.read_param(self.parameters)
        
        return