"""
This is an example of a post-processing script. This will read in the input and output files and process them like as in
a typical CTEM simulation

"""
import ctem
sim = ctem.Simulation(isnew=False)
sim.read_output()
sim.user['ncount'] += 1 # Increment the output count in the situation in which the simulation was run outside of the driver script
sim.process_output()