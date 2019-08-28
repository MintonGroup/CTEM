import numpy as np
import pandas as pd

gridsize = 2000
pix = 3.08e3
side = (gridsize * pix)
area = side**2
iside = [-side, 0.0, side]
basin_crater_cutoff = 300000.0 #diameter in meters

craters = pd.read_csv('NPF-global-Kd0.0001/dist/ocum_000100.dat',delim_whitespace=True)

basins = craters.loc[craters['#Dcrat(m)'] >= basin_crater_cutoff]

for bind, b in basins.iterrows():
    print(f"Basin{bind:02}: ",b['#Dcrat(m)'],b['time(y)'])
    overlaps = open(f"NPF-global-Kd0.0001/basin_overlap/overlap{bind:02}.dat", 'w')
    print(f"#BASIN NUMBER {bind}: Diameter = {b['#Dcrat(m)']}, Time {b['time(y)']}", file=overlaps)
    for cind, c in craters.iterrows():
        if c['time(y)'] > b['time(y)']:
            for i in iside:
                for j in iside:
                    disx = b['xpos(m)'] - c['xpos(m)'] + i
                    disy = b['ypos(m)'] - c['ypos(m)'] + j
                    distance = np.sqrt(disx ** 2 + disy ** 2)
                    if distance < 0.5 * (b['#Dcrat(m)'] + c['#Dcrat(m)']):
                        print(c['#Dcrat(m)'],c['xpos(m)'],c['ypos(m)'],c['time(y)'], file=overlaps)


