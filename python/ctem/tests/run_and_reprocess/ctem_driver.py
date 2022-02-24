"""
This is an example of a standard "driver" script. This will read in the input files, pre-process them for the main
Fortran code and post-process them.
"""
import ctem
sim = ctem.Simulation()
sim.run()