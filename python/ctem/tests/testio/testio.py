import ctem

d = dir(ctem)

sim = ctem.Simulation()
for k, v in sim.parameters.items():
    print(k, v)