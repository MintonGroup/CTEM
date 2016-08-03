#!/usr/local/bin/python
#
#Cratered Terrain Evolution Model file reading utilities
#
#Original IDL design: Jim Richardson, Arecibo Observatory
#Revised IDL design: David Minton, Purdue University
#Re-engineered design and Python implementation: Matthew Route, Purdue University
#August 2016
#
#Known issues for operation with CTEM
#1) ctem.in has 1.0d0 value for maxcrat which is not readable by Python

import numpy

def read_ctemin(parameters,notset):
    #Read and parse ctem.in file
   inputfile = parameters['workingdir'] + parameters['ctemfile']

   #Read ctem.in file
   print 'Reading input file '+ parameters['ctemfile']
   fp = open(inputfile,'r')
   lines = fp.readlines()
   fp.close()
   
   #Process file text
   for line in lines:
      fields = line.split()
      if len(fields) > 0:
         if ('pix' == fields[0].lower()): parameters['pix']=float(fields[1])
         if ('gridsize' == fields[0].lower()): parameters['gridsize']=int(fields[1])
         if ('seed' == fields[0].lower()): parameters['seed']=int(fields[1])
         if ('sfdfile' == fields[0].lower()): parameters['sfdfile']=fields[1]
         if ('impfile' == fields[0].lower()): parameters['impfile']=fields[1]
         if ('maxcrat' == fields[0].lower()): parameters['maxcrat']=float(fields[1])
         if ('sfdcompare' == fields[0].lower()): parameters['sfdcompare']=fields[1]
         if ('interval' == fields[0].lower()): parameters['interval']=float(fields[1])
         if ('numintervals' == fields[0].lower()): parameters['numintervals']=int(fields[1])
         if ('popupconsole' == fields[0].lower()): parameters['popupconsole']=fields[1]
         if ('saveshaded' == fields[0].lower()): parameters['saveshaded']=fields[1]
         if ('saverego' == fields[0].lower()): parameters['saverego']=fields[1]
         if ('savepres' == fields[0].lower()): parameters['savepres']=fields[1]
         if ('savetruelist' == fields[0].lower()): parameters['savetruelist']=fields[1]
         if ('runtype' == fields[0].lower()): parameters['runtype']=fields[1]
         if ('restart' == fields[0].lower()): parameters['restart']=fields[1]
         if ('shadedminh' == fields[0].lower()):
            parameters['shadedminh'] = float(fields[1])
            parameters['shadedminhdefault'] = 0
         if ('shadedmaxh' == fields[0].lower()):
            parameters['shadedmaxh'] = float(fields[1])
            parameters['shadedmaxhdefault'] = 0
       
   #Test values for further processing
   if (parameters['interval'] <= 0.0):
       print 'Invalid value for or missing variable INTERVAL in '+ inputfile
   if (parameters['numintervals'] <= 0):
       print 'Invalid value for or missing variable NUMINTERVALS in '+ inputfile
   if (parameters['pix'] <= 0.0):
       print 'Invalid value for or missing variable PIX in '+ inputfile
   if (parameters['gridsize'] <= 0):
       print 'Invalid value for or missing variable GRIDSIZE in '+ inputfile
   if (parameters['seed'] == 0):
       print 'Invalid value for or missing variable SEED in '+ inputfile
   if (parameters['sfdfile'] == notset):
       print 'Invalid value for or missing variable SFDFILE in '+ inputfile
   if (parameters['impfile'] == notset):
       print 'Invalid value for or missing variable IMPFILE in '+ inputfile
   if (parameters['popupconsole'] == notset):
       print 'Invalid value for or missing variable POPUPCONSOLE in '+ inputfile
   if (parameters['saveshaded'] == notset):
       print 'Invalid value for or missing variable SAVESHADED in '+ inputfile
   if (parameters['saverego'] == notset):
       print 'Invalid value for or missing variable SAVEREGO in '+ inputfile
   if (parameters['savepres'] == notset):
       print 'Invalid value for or missing variable SAVEPRES in '+ inputfile
   if (parameters['savetruelist'] == notset):
       print 'Invalid value for or missing variable SAVETRUELIST in '+ inputfile
   if (parameters['runtype'] == notset):
       print 'Invalid value for or missing variable RUNTYPE in '+ inputfile
   if (parameters['restart'] == notset):
       print 'Invalid value for or missing variable RESTART in '+ inputfile
       
   return

def read_formatted_ascii(filename, skip_lines):
   #Generalized ascii text reader
   #For use with sfdcompare, production, odist, tdist, pdist data files
   data = numpy.genfromtxt(filename, skip_header = skip_lines)
   return data

def read_unformatted_binary(filename, gridsize):
   #Read unformatted binary files created by Fortran
   #For use with surface ejecta and surface dem data files
   dt = numpy.float
   data = numpy.fromfile(filename, dtype = dt)
   data.shape = (gridsize,gridsize)
   
   return data

def read_ctemdat(parameters, seedarr):
    #Read and parse ctem.dat file
    datfile = parameters['workingdir'] + 'ctem.dat'
 
    #Read ctem.dat file    
    print 'Reading input file '+ parameters['datfile']
    fp = open(datfile,'r')
    lines = fp.readlines()
    fp.close()
   
    #Parse file lines and update parameter fields
    fields = lines[0].split()
    if len(fields) > 0:
      parameters['totalimpacts'] = float(fields[0])
      parameters['ncount'] = int(fields[1])
      parameters['curyear'] = float(fields[2])
      parameters['restart'] = fields[3]
      parameters['fracdone'] = float(fields[4])
      parameters['masstot'] = float(fields[5])

    #Parse remainder of file to build seed array
    nlines = len(lines)
    index = 1
    while (index < nlines):
        fields = lines[index].split()
        seedarr[index - 1] = float(fields[0])
        index += 1

    parameters['seedn'] = index - 1

    return

def read_impact_mass(filename):
   #Read impact mass file
   
   fp=open(filename,'r')
   line=fp.readlines()
   fp.close()
   
   fields = line[0].split()
   if (len(fields) > 0):
      mass = float(fields[0])
   else:
      mass = 0
      
   return mass