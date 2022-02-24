"""
This is an example of a post-processing script. This will read in the input and output files and process them like as in
a typical CTEM simulation

"""
import ctem
sim = ctem.Simulation(isnew=False)
sim.read_output()
sim.user['ncount'] += 1
sim.process_output()