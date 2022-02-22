import numpy as np
import os
import subprocess
import shutil
from ctem import io

class Simulation:
    """
    This is a class that defines the basic CTEM simulation object. It initializes a dictionary of user and reads
    it in from a file (default name is ctem.in). It also creates the directory structure needed to store simulation
    files if necessary.
    """
    def __init__(self, param_file="ctem.in"):
        currentdir = os.getcwd()
        self.user = {
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
            'ctemfile': os.path.join(currentdir, param_file),
            'impfile': None,
            'sfdcompare': None,
            'sfdfile': None
        }
        self.user = io.read_user_input(self.user)
        
        self.output_filenames = {
            'dat': 'ctem.dat',
            'dem': 'surface_dem.dat',
            'diam': 'surface_diam.dat',
            'ejc': 'surface_ejc.dat',
            'pos': 'surface_pos.dat',
            'time': 'surface_time.dat',
            'ocum': 'ocumulative.dat',
            'odist': 'odistribution.dat',
            'pdist': 'pdistribution.dat',
            'tcum' : 'tcumulative.dat',
            'tdist' : 'tdistribution.dat',
            'impmass' : 'impactmass.dat',
            'fracdone' : 'fracdone.dat',
            'regodepth' : 'regolithdepth.dat',
        }
        
        for k, v in self.output_filenames.items():
            self.output_filenames[k] = os.path.join(currentdir, v)
        
        # Set up data arrays
        self.seedarr = np.zeros(100, dtype=int)
        self.seedarr[0] = self.user['seed']
        self.odist = np.zeros([1, 6])
        self.pdist = np.zeros([1, 6])
        self.tdist = np.zeros([1, 6])
        self.surface_dem = np.zeros([self.user['gridsize'], self.user['gridsize']], dtype=float)
        self.regolith = np.zeros([self.user['gridsize'], self.user['gridsize']], dtype=float)
        
        if self.user['sfdcompare'] is not None:
            # Read sfdcompare file
            sfdfile = os.path.join(self.user['workingdir'], self.user['sfdcompare'])
            self.ph1 = io.read_formatted_ascii(sfdfile, skip_lines=0)
           
        # Scale the production function to the simulation domain
        self.scale_production()

        # Starting new or old run?
        if (self.user['restart'].upper() == 'F'):
            print('Starting a new run')

            io.create_dir_structure(self.user)
            
            if (self.user['runtype'].upper() == 'STATISTICAL'):
                self.user['ncount'] = 1

                # Write ctem.dat file
                io.write_datfile(self.user, self.output_filenames['dat'], self.seedarr)
            else:
                self.user['ncount'] = 0

            # Delete any old output files
            for k, v in self.output_filenames.items():
                if os.path.isfile(v):
                    os.remove(v)
        else:
            print('Continuing a previous run')
            self.process_interval(isnew=False)


        return
    
    def scale_production(self):
        """
        Scales the production function to the simulation domain and saves the scaled production function to a file that
        will be read in by the main Fortran program.
        """
        
        # Read production function file
        impfile = os.path.join(self.user['workingdir'], self.user['impfile'])
        prodfunction = io.read_formatted_ascii(impfile, skip_lines=0)

        # Create impactor production population
        area = (self.user['gridsize'] * self.user['pix']) ** 2
        self.production = np.copy(prodfunction)
        self.production[:, 1] = self.production[:, 1] * area * self.user['interval']

        # Write the scaled production function to file
        io.write_production(self.user, self.production)
        
        return
    
    def run(self):
        """
        Runs a complete simulation over all intervals specified in the input user. This method loops over all
        intervals and process outputs into images, and redirect output files to their apporpriate folders.
        
        This method replaces most of the functionality of the original ctem_driver.py
        """

        print('Beginning loops')
        while (self.user['ncount'] <= self.user['numintervals']):
            if (self.user['ncount'] > 0):
                # Move ctem.dat
                self.redirect_outputs(['dat'], 'misc')
                
                #Execute Fortran code
                self.compute_one_interval()
                
                # Process output files
                self.process_interval()

            # Display results
            print(self.user['ncount'], '  Displaying results')

            # Write surface dem, surface ejecta, shaded relief, and rplot data
            io.image_dem(self.user, self.surface_dem)
            if (self.user['saverego'].upper() == 'T'):
                io.image_regolith(self.user, self.surface_ejc)
            if (self.user['saveshaded'].upper() == 'T'):
                io.image_shaded_relief(self.user, self.surface_dem)
            if (self.user['savepres'].upper() == 'T'):
                io.create_rplot(self.user, self.odist, self.pdist, self.tdist, self.ph1)

            # Update ncount
            self.user['ncount'] = self.user['ncount'] + 1

            if ((self.user['runtype'].upper() == 'STATISTICAL') or (self.user['ncount'] == 1)):
                self.user['restart'] = 'F'
                self.user['curyear'] = 0.0
                self.user['totalimpacts'] = 0
                self.user['masstot'] = 0.0

                # Delete tdistribution file, if it exists
                tdist_file = self.user['workingdir'] + 'tdistribution.dat'
                if os.path.isfile(tdist_file):
                    os.remove(tdist_file)

            else:
                self.user['restart'] = 'T'

            # Write ctem.dat file
            io.write_datfile(self.user, self.output_filenames['dat'],self.seedarr)

        return
    
    def compute_one_interval(self):
        """
        Executes the Fortran code to generate one interval of simulation output.
        """
        # Create crater population and display CTEM progress on screen
        print(self.user['ncount'], '  Calling FORTRAN routine')
        with subprocess.Popen([os.path.join(self.user['workingdir'], 'CTEM')],
                              stdout=subprocess.PIPE,
                              universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')

        return
    
    def process_interval(self, isnew=True):
        """
        Reads in all of the files generated by the Fortran code after one interval and redirects them to appropriate
        folders for storage after appending the interval number to the filename.
        """
        # Read Fortran output
        print(self.user['ncount'], '  Reading Fortran output')

        # Read surface dem(shaded relief) and ejecta data files
        self.surface_dem = io.read_unformatted_binary(self.output_filenames['dem'], self.user['gridsize'])
        self.surface_ejc = io.read_unformatted_binary(self.output_filenames['ejc'], self.user['gridsize'])

        # Read odistribution, tdistribution, and pdistribution files
        self.odist = io.read_formatted_ascii(self.output_filenames['odist'], skip_lines=1)
        self.tdist = io.read_formatted_ascii(self.output_filenames['tdist'], skip_lines=1)
        self.pdist = io.read_formatted_ascii(self.output_filenames['pdist'], skip_lines=1)

        # Read impact mass from file
        impact_mass = io.read_impact_mass(self.output_filenames['impmass'])

        # Read ctem.dat file
        io.read_datfile(self.user, self.output_filenames['dat'], self.seedarr)

        # Save copy of crater distribution files
        # Update user: mass, curyear, regolith properties
        self.user['masstot'] = self.user['masstot'] + impact_mass

        self.user['curyear'] = self.user['curyear'] + self.user['fracdone'] * self.user['interval']
        template = "%(fracdone)9.6f %(curyear)19.12E\n"
        with open(self.output_filenames['fracdone'], 'w') as fp_frac:
            fp_frac.write(template % self.user)

        reg_text = "%19.12E %19.12E %19.12E %19.12E\n" % (self.user['curyear'],
                                                          np.mean(self.surface_ejc), np.amax(self.surface_ejc),
                                                          np.amin(self.surface_ejc))
        with open(self.output_filenames['regodepth'], 'w') as fp_reg:
            fp_reg.write(reg_text)

        if isnew:
            self.redirect_outputs(['odist', 'ocum', 'pdist', 'tdist'], 'dist')
            if (self.user['savetruelist'].upper() == 'T'):
                self.redirect_outputs(['tcum'], 'dist')
            self.redirect_outputs(['impmass'], 'misc')
        return

    def redirect_outputs(self, filekeys, foldername):
        """
        Copies a set of output files from the working directory into a subfolder.
        Takes as input the dictionaries of user parameters and output file names, the dictionary keys to the file names
        you wish to redirect, and the redirection destination folder name. The new name will be the key name + zero padded
        interval number.

        Example:
            calling redirect_outputs(['odist','tdist'], 'dist') when user['ncount'] is 1
            should do the following:
            copy 'odistribution.dat' to 'dist/odist_000001.dat'
            copy 'tdistribution.dat' to 'dist/tdist_000001.dat'
        """
        
        for k in filekeys:
            forig = self.output_filenames[k]
        fdist = os.path.join(self.user['workingdir'], foldername, f"{k}_{self.user['ncount']:06d}.dat")

if __name__ == '__main__':
    sim = Simulation()
    sim.run()