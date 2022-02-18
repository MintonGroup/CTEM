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
ctem_io_writers.create_dir_structure(parameters)

dem_file = parameters['workingdir'] + 'surface_dem.dat'
surface_dem = ctem_io_readers.read_unformatted_binary(dem_file, parameters['gridsize'])

#Set up data arrays
seedarr = numpy.zeros(100, dtype = numpy.int)
seedarr[0] = parameters['seed']

#Read ctem.dat file
ctem_io_readers.read_ctemdat(parameters, seedarr)

#Write surface dem
ctem_io_writers.image_dem(parameters, surface_dem)