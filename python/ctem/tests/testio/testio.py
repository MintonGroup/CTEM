import ctem

d = dir(ctem)

sim = ctem.Simulation()
for k, v in sim.user.items():
    print(k, v)