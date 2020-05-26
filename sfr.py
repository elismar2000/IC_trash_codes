import matplotlib.pyplot as plt
import numpy as np
from pygadgetreader import *

snap0 = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test1/snapshot_0000'
snap1 = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test1/snapshot_0001'

pos_gas0 = readsnap(snap0, 'pos', 'gas')
pos_gas1 = readsnap(snap1, 'pos', 'gas')
pos_star1 = readsnap(snap1, 'pos', 'star')

id_gas0 = readsnap(snap0, 'pid', 'gas')
id_gas1 = readsnap(snap1, 'pid', 'gas')
id_star1 = readsnap(snap1, 'pid', 'star')

print('pos_gas0: ', len(pos_gas0[:, 0]))
print('pos_gas1: ', len(pos_gas1[:, 0]))
print('pos_star1: ', len(pos_star1[:, 0]))

missing_gas = np.isin(id_gas0, id_gas1)

plt.plot(pos_gas0[:, 0], pos_gas0[:, 1], 'r.')
plt.plot(pos_gas1[:, 0], pos_gas1[:, 1], 'b.')
plt.plot(pos_star1[:, 0], pos_star1[:, 1], 'y*')
plt.show()
