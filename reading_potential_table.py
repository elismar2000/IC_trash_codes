from tables import *
import numpy as np
import matplotlib.pyplot as plt

h5file = open_file('galmer-like_sim_minima.h5', 'a')
table = h5file.root.potential.readout

xmin1 = np.array([j['xmin1'] for j in table.iterrows()])
ymin1 = np.array([j['ymin1'] for j in table.iterrows()])

xmin2 = np.array([j['xmin2'] for j in table.iterrows()])
ymin2 = np.array([j['ymin2'] for j in table.iterrows()])

snapshots = np.array([ j['snapshot'] for j in table.iterrows()])
print(snapshots)

plt.plot(xmin1, ymin1, '-o', label='gal1')
plt.plot(xmin2, ymin2, '-o', label='gal2')
plt.legend()
plt.show()
