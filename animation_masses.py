import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import Animation
from potential_energy import select_particles

fig, ax = plt.subplots(figsize = (8, 8))

table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"

snapshot = 20
galaxy = 1

t = Table.read(table_path, snapshot)
mask = select_particles(table_path, snapshot, galaxy)
x = t['X'][mask]
y = t['Y'][mask]
z = t['Z'][mask]
coords = np.transpose(np.array([x, y, z]))
masses, edges = np.histogramdd(coords, bins = 25, weights = t['MASS'][mask])

im = plt.imshow(masses[0])

def animate(bin):
    t = Table.read(table_path, snapshot)
    mask = select_particles(table_path, snapshot, galaxy)
    x = t['X'][mask]
    y = t['Y'][mask]
    z = t['Z'][mask]
    coords = np.transpose(np.array([x, y, z]))
    masses, edges = np.histogramdd(coords, bins = 25, weights = t['MASS'][mask])
    im.set_data(masses[bin])

ani = animation.FuncAnimation(fig, animate, np.arange(0, 26, 1), interval = 100, blit = False)

plt.show()
