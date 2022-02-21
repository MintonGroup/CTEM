import numpy as np
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
        
        # Set up data arrays
        self.seedarr = np.zeros(100, dtype=int)
        self.seedarr[0] = self.parameters['seed']
        self.odist = np.zeros([1, 6])
        self.pdist = np.zeros([1, 6])
        self.tdist = np.zeros([1, 6])
        self.surface_dem = np.zeros([self.parameters['gridsize'], self.parameters['gridsize']], dtype=np.float)
        self.regolith = np.zeros([self.parameters['gridsize'], self.parameters['gridsize']], dtype=np.float)
        
        if self.parameters['sfdcompare'] is not None:
            # Read sfdcompare file
            sfdfile = os.path.join(self.parameters['workingdir'], self.parameters['sfdcompare'])
            self.ph1 = io.read_formatted_ascii(sfdfile, skip_lines=0)

        # Read production function file
        impfile = os.path.join(self.parameters['workingdir'], parameters['impfile'])
        prodfunction = io.read_formatted_ascii(impfile, skip_lines=0)

        # Create impactor production population
        area = (self.parameters['gridsize'] * self.parameters['pix']) ** 2
        self.production = np.copy(prodfunction)
        self.production[:, 1] = self.production[:, 1] * area * self.parameters['interval']

        # Write corrected production function to file
        io.write_production(self.parameters, self.production)

        # Starting new or old run?
        if (self.parameters['restart'].upper() == 'F'):
            print('Starting a new run')
    
            if (self.parameters['runtype'].upper() == 'STATISTICAL'):
                self.parameters['ncount'] = 1

                # Write ctem.dat file
                io.write_ctemdat(self.parameters, self.seedarr)
    
            else:
                self.parameters['ncount'] = 0

            # Delete tdistribution file, if it exists
            tdist_file = os.path.join(self.parameters['workingdir'], 'tdistribution.dat')
            if os.path.isfile(tdist_file):
                os.remove(tdist_file)
        else:
            print('Continuing a previous run')

            # Read surface dem(shaded relief) and ejecta data files
            dem_file = os.path.join(self.parameters['workingdir'], 'surface_dem.dat')
            self.surface_dem = io.read_unformatted_binary(dem_file, self.parameters['gridsize'])
            ejecta_file = os.path.join(self.parameters['workingdir'], 'surface_ejc.dat')
            self.regolith = io.read_unformatted_binary(ejecta_file, self.parameters['gridsize'])

            # Read odistribution, tdistribution, and pdistribution files
            ofile = os.path.join(self.parameters['workingdir'],'odistribution.dat')
            odist = io.read_formatted_ascii(ofile, skip_lines=1)
            tfile = os.path.join(self.parameters['workingdir'], 'tdistribution.dat')
            tdist = io.read_formatted_ascii(tfile, skip_lines=1)
            pfile = os.path.join(self.parameters['workingdir'], 'pdistribution.dat')
            pdist = io.read_formatted_ascii(pfile, skip_lines=1)

            # Read impact mass from file
            massfile = os.path.join(self.parameters['workingdir'], 'impactmass.dat')
            impact_mass = io.read_impact_mass(massfile)

            # Read ctem.dat file
            io.read_ctemdat(self.parameters, self.seedarr)

            # Open fracdonefile and regodepthfile for writing
        # filename = os.path.join(self.parameters['workingdir'], 'fracdone.dat')
        # fp_frac = open(filename, 'w')
        # filename = os.path.join(self.parameters['workingdir'], 'regolithdepth.dat')
        # fp_reg = open(filename, 'w')

        io.create_dir_structure(self.parameters)


        return
    
    
    def