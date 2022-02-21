import numpy
import os

def real2float(realstr):
    """
    Converts a Fortran-generated ASCII string of a real value into a numpy float type. Handles cases where double precision
    numbers in exponential notation use 'd' or 'D' instead of 'e' or 'E'

    Parameters
    ----------
    realstr : string
        Fortran-generated ASCII string of a real value.

    Returns
    -------
     : float
        The converted floating point value of the input string
    """
    return float(realstr.replace('d', 'E').replace('D', 'E'))

def read_param(parameters):
    # Read and parse ctem.in file
    inputfile = os.path.join(parameters['workingdir'], parameters['ctemfile'])

    # Read ctem.in file
    print('Reading input file ' + parameters['ctemfile'])
    with open(inputfile, 'r') as fp:
        lines = fp.readlines()

    # Process file text
    for line in lines:
        fields = line.split()
        if len(fields) > 0:
            if ('pix' == fields[0].lower()): parameters['pix'] = real2float(fields[1])
            if ('gridsize' == fields[0].lower()): parameters['gridsize'] = int(fields[1])
            if ('seed' == fields[0].lower()): parameters['seed'] = int(fields[1])
            if ('sfdfile' == fields[0].lower()): parameters['sfdfile'] = fields[1]
            if ('impfile' == fields[0].lower()): parameters['impfile'] = fields[1]
            if ('maxcrat' == fields[0].lower()): parameters['maxcrat'] = real2float(fields[1])
            if ('sfdcompare' == fields[0].lower()): parameters['sfdcompare'] = fields[1]
            if ('interval' == fields[0].lower()): parameters['interval'] = real2float(fields[1])
            if ('numintervals' == fields[0].lower()): parameters['numintervals'] = int(fields[1])
            if ('popupconsole' == fields[0].lower()): parameters['popupconsole'] = fields[1]
            if ('saveshaded' == fields[0].lower()): parameters['saveshaded'] = fields[1]
            if ('saverego' == fields[0].lower()): parameters['saverego'] = fields[1]
            if ('savepres' == fields[0].lower()): parameters['savepres'] = fields[1]
            if ('savetruelist' == fields[0].lower()): parameters['savetruelist'] = fields[1]
            if ('runtype' == fields[0].lower()): parameters['runtype'] = fields[1]
            if ('restart' == fields[0].lower()): parameters['restart'] = fields[1]
            if ('shadedminh' == fields[0].lower()):
                parameters['shadedminh'] = real2float(fields[1])
                parameters['shadedminhdefault'] = 0
            if ('shadedmaxh' == fields[0].lower()):
                parameters['shadedmaxh'] = real2float(fields[1])
                parameters['shadedmaxhdefault'] = 0

    # Test values for further processing
    if (parameters['interval'] <= 0.0):
        print('Invalid value for or missing variable INTERVAL in ' + inputfile)
    if (parameters['numintervals'] <= 0):
        print('Invalid value for or missing variable NUMINTERVALS in ' + inputfile)
    if (parameters['pix'] <= 0.0):
        print('Invalid value for or missing variable PIX in ' + inputfile)
    if (parameters['gridsize'] <= 0):
        print('Invalid value for or missing variable GRIDSIZE in ' + inputfile)
    if (parameters['seed'] == 0):
        print('Invalid value for or missing variable SEED in ' + inputfile)
    if (parameters['sfdfile'] is None):
        print('Invalid value for or missing variable SFDFILE in ' + inputfile)
    if (parameters['impfile'] is None):
        print('Invalid value for or missing variable IMPFILE in ' + inputfile)
    if (parameters['popupconsole'] is None):
        print('Invalid value for or missing variable POPUPCONSOLE in ' + inputfile)
    if (parameters['saveshaded'] is None):
        print('Invalid value for or missing variable SAVESHADED in ' + inputfile)
    if (parameters['saverego'] is None):
        print('Invalid value for or missing variable SAVEREGO in ' + inputfile)
    if (parameters['savepres'] is None):
        print('Invalid value for or missing variable SAVEPRES in ' + inputfile)
    if (parameters['savetruelist'] is None):
        print('Invalid value for or missing variable SAVETRUELIST in ' + inputfile)
    if (parameters['runtype'] is None):
        print('Invalid value for or missing variable RUNTYPE in ' + inputfile)
    if (parameters['restart'] is None):
        print('Invalid value for or missing variable RESTART in ' + inputfile)
    
    return parameters
