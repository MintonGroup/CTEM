#!/usr/local/bin/python
#
#Cratered Terrain Evolution Model driver
#
#Original IDL design: Jim Richardson, Arecibo Observatory
#Revised IDL design: David Minton, Purdue University
#Re-engineered design and Python implementation: Matthew Route, Purdue University
#August 2016

#Import general purpose modules
import numpy
import os
import subprocess
import shutil

#Import CTEM modules
import ctem_io_readers
import ctem_io_writers

#Create and initialize data dictionaries for parameters and options from CTEM.in
notset = '-NOTSET-'
currentdir = os.getcwd() + os.sep

parameters={'restart': notset,
            'runtype': notset,
            'popupconsole': notset,
            'saveshaded': notset,
            'saverego': notset,
            'savepres': notset,
            'savetruelist': notset,
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
            'ctemfile': 'ctem.in',
            'datfile': 'ctem.dat',
            'impfile': notset,
            'sfdcompare': notset,
            'sfdfile': notset}

#Read ctem.in to initialize parameter values based on user input
ctem_io_readers.read_ctemin(parameters,notset)

#Read sfdcompare file
sfdfile = parameters['workingdir'] + parameters['sfdcompare']
ph1 = ctem_io_readers.read_formatted_ascii(sfdfile, skip_lines = 0)

#Set up data arrays
seedarr = numpy.zeros(100, dtype = numpy.int)
seedarr[0] = parameters['seed']
odist = numpy.zeros([1, 6])
pdist = numpy.zeros([1, 6])
tdist = numpy.zeros([1, 6])
surface_dem = numpy.zeros([parameters['gridsize'], parameters['gridsize']], dtype = numpy.float)
regolith = numpy.zeros([parameters['gridsize'], parameters['gridsize']], dtype =numpy.float)

#Read production function file
impfile = parameters['workingdir'] + parameters['impfile']
prodfunction = ctem_io_readers.read_formatted_ascii(impfile, skip_lines = 0)

#Create impactor production population
area = (parameters['gridsize'] * parameters['pix'])**2
production = numpy.copy(prodfunction)
production[:,1] = production[:,1] * area * parameters['interval']

#Write corrected production function to file
ctem_io_writers.write_production(parameters, production)

#Starting new or old run?
if (parameters['restart'].upper() == 'F'):
    print('Starting a new run')
    
    if (parameters['runtype'].upper() == 'STATISTICAL'):
        parameters['ncount'] = 1
        
        #Write ctem.dat file
        ctem_io_writers.write_ctemdat(parameters, seedarr)
        
    else:
        parameters['ncount'] = 0
    
    #Delete tdistribution file, if it exists
    tdist_file = parameters['workingdir'] + 'tdistribution.dat'
    if os.path.isfile(tdist_file):
        os.remove(tdist_file)
    
else:
    print('Continuing a previous run')
    
    #Read surface dem(shaded relief) and ejecta data files
    dem_file = parameters['workingdir'] + 'surface_dem.dat'
    surface_dem = ctem_io_readers.read_unformatted_binary(dem_file, parameters['gridsize'])
    ejecta_file = parameters['workingdir'] + 'surface_ejc.dat'
    regolith = ctem_io_readers.read_unformatted_binary(ejecta_file, parameters['gridsize'])
    
    #Read odistribution, tdistribution, and pdistribution files
    ofile = parameters['workingdir'] + 'odistribution.dat'
    odist = ctem_io_readers.read_formatted_ascii(ofile, skip_lines = 1)
    tfile = parameters['workingdir'] + 'tdistribution.dat'
    tdist = ctem_io_readers.read_formatted_ascii(tfile, skip_lines = 1)
    pfile = parameters['workingdir'] + 'pdistribution.dat'
    pdist = ctem_io_readers.read_formatted_ascii(pfile, skip_lines = 1)
    
    #Read impact mass from file
    massfile = parameters['workingdir'] + 'impactmass.dat'
    impact_mass = ctem_io_readers.read_impact_mass(massfile)
    
    #Read ctem.dat file
    ctem_io_readers.read_ctemdat(parameters, seedarr)  
    
#Open fracdonefile and regodepthfile for writing
filename = parameters['workingdir'] + 'fracdone.dat'
fp_frac = open(filename,'w')
filename = parameters['workingdir'] + 'regolithdepth.dat'
fp_reg = open(filename,'w')

#Begin CTEM processing loops
print('Beginning loops')

ctem_io_writers.create_dir_structure(parameters)

while (parameters['ncount'] <= parameters['numintervals']):
    
    #Create crater population
    if (parameters['ncount'] > 0):
        
        #Move ctem.dat
        forig = parameters['workingdir'] + 'ctem.dat'
        fdest = parameters['workingdir'] + 'misc' + os.sep + "ctem%06d.dat" % parameters['ncount']
        shutil.copy2(forig, fdest)
    
        #Create crater population and display CTEM progress on screen
        print(parameters['ncount'], '  Calling FORTRAN routine')
        with subprocess.Popen([parameters['workingdir']+'CTEM'], stdout=subprocess.PIPE, bufsize=1,universal_newlines=True) as p:
            for line in p.stdout:
                print(line, end='')

        
        #Optional: do not pipe CTEM progress to the screen
        #subprocess.check_output([parameters['workingdir']+'CTEM'])
    
        #Read Fortran output
        print(parameters['ncount'], '  Reading Fortran output')
    
        #Read surface dem(shaded relief) and ejecta data files
        dem_file = parameters['workingdir'] + 'surface_dem.dat'
        surface_dem = ctem_io_readers.read_unformatted_binary(dem_file, parameters['gridsize'])    
        ejecta_file = parameters['workingdir'] + 'surface_ejc.dat'
        regolith = ctem_io_readers.read_unformatted_binary(ejecta_file, parameters['gridsize'])
    
        #Read odistribution, tdistribution, and pdistribution files
        ofile = parameters['workingdir'] + 'odistribution.dat'
        odist = ctem_io_readers.read_formatted_ascii(ofile, skip_lines = 1)
        tfile = parameters['workingdir'] + 'tdistribution.dat'
        tdist = ctem_io_readers.read_formatted_ascii(tfile, skip_lines = 1)
        pfile = parameters['workingdir'] + 'pdistribution.dat'
        pdist = ctem_io_readers.read_formatted_ascii(pfile, skip_lines = 1)
    
        #Read impact mass from file
        massfile = parameters['workingdir'] + 'impactmass.dat'
        impact_mass = ctem_io_readers.read_impact_mass(massfile)
    
        #Read ctem.dat file
        ctem_io_readers.read_ctemdat(parameters, seedarr)
        
        #Update parameters: mass, curyear, regolith properties
        parameters['masstot'] = parameters['masstot'] + impact_mass
        
        parameters['curyear'] = parameters['curyear'] + parameters['fracdone'] * parameters['interval']
        template = "%(fracdone)9.6f %(curyear)19.12E\n"
        fp_frac.write(template % parameters)
        
        reg_text = "%19.12E %19.12E %19.12E %19.12E\n" % (parameters['curyear'],
                numpy.mean(regolith), numpy.amax(regolith), numpy.amin(regolith))
        fp_reg.write(reg_text)
        
        #Save copy of crater distribution files
        ctem_io_writers.copy_dists(parameters)

    #Display results
    print(parameters['ncount'], '  Displaying results')

    #Write surface dem, surface ejecta, shaded relief, and rplot data
    ctem_io_writers.image_dem(parameters, surface_dem)
    if (parameters['saverego'].upper() == 'T'):
        ctem_io_writers.image_regolith(parameters, regolith)
    if (parameters['saveshaded'].upper() == 'T'):
        ctem_io_writers.image_shaded_relief(parameters, surface_dem)
    if (parameters['savepres'].upper() == 'T'):
        ctem_io_writers.create_rplot(parameters,odist,pdist,tdist,ph1)
     
    #Update ncount
    parameters['ncount'] = parameters['ncount'] + 1

    if ((parameters['runtype'].upper() == 'STATISTICAL') or (parameters['ncount'] == 1)):
        parameters['restart'] = 'F'
        parameters['curyear'] = 0.0
        parameters['totalimpacts'] = 0
        parameters['masstot'] = 0.0
        
        #Delete tdistribution file, if it exists
        tdist_file = parameters['workingdir'] + 'tdistribution.dat'
        if os.path.isfile(tdist_file):
            os.remove(tdist_file)
            
    else:
        parameters['restart'] = 'T'
    
    #Write ctem.dat file
    ctem_io_writers.write_ctemdat(parameters, seedarr)    

#Close updateable fracdonefile and regodepthfile files
fp_frac.close()
fp_reg.close()