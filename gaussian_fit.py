from scipy.optimize import curve_fit, minimize
import numpy as np
import matplotlib.pyplot as plt
from potential_3d import Potential

p = Potential()

t_path = "/home/elismar/Documentos/Fisica/IC/GalMer/inflow/tables_arp142_v2"
snapshot = 71

coords_min, pot, c, i = p.pot3d(t_path, snapshot, n_bins=20, smooth=1.2, width=1.5)
pot_x = pot[:, int(i[1]), int(i[2])]

def func(x, mu, sigma):
    return (1/(np.sqrt(2*np.pi*sigma**2)))*np.exp(-(x-mu)**2/(2*sigma**2))

popt, popv = curve_fit(func, c[0], pot_x)

g = func(c[0], popt[0], popt[1])

plt.plot(c[0], pot_x, '--', label='Discrete potential')
plt.plot(c[0], g, label='best fit')
