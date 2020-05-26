import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import Animation
from simulation_viewer import bins_mask_again

fig, ax = plt.subplots(1, 3, figsize = (10, 10), tight_layout=True)
#fig = plt.figure()

table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"

t = Table.read(table_path, 1)
data, xedges, yedges = np.histogram2d(t['X'], t['Y'], bins=200)
im = plt.imshow(data)

def animate(snapshot):
    t = Table.read(table_path, snapshot)
    x = t['X']
    y = t['Y']
    z = t['Z']
    bins = 200
    fig1 = ax[0].hist2d(x, y, bins=bins, weights=t['MASS'])
    fig2 = ax[1].hist2d(x, z, bins=bins, weights=t['MASS'])
    fig3 = ax[2].hist2d(y, z, bins=bins, weights=t['MASS'])
    fig.tight_layout()

ani = animation.FuncAnimation(fig, animate, np.arange(1, 72, 1), interval = 100, blit = False)

#ani.save('x_vs_y(1).mp4')

plt.show()
