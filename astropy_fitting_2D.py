from potential_2d import Potential
import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from inflow import Inflow
from astropy.table import Table

p = Potential()
t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/tables_arp142_v2'
snapshot = 1
pot, c, d_bin = p.pot2d(t_path, snapshot, zoom=0)

g_init = models.Gaussian1D(amplitude=-1., mean=0, stddev=0.5)
fit_g = fitting.LevMarLSQFitter()

def evaluate_minimum(dimension: str = 'x'):

    if dimension == 'x':
        mins = np.zeros_like(c[0])
        values = np.zeros_like(c[0])
        for i in range(np.size(pot, axis=1)):
            g = fit_g(g_init, c[0], pot[:, i])
            m = minimize(g, 0)
            mins[i] = m.x
            values[i] = m.fun

        i1 = np.where(values == values.min())
        values = np.delete(values, i1)
        i2 = np.where(values == values.min())

        x1 = mins_x[i1]
        x2 = mins_x[i2]
        y1 = c[1, i1]
        y2 = c[1, i2]

        coords_min1 = np.array([x1, y1])
        coords_min2 = np.array([x2, y2])

    if dimension == 'y':
        mins = np.zeros_like(c[1])
        values = np.zeros_like(c[1])
        for i in range(np.size(pot, axis=0)):
            g = fit_g(g_init, c[1], pot[i])
            m = minimize(g, 0)
            mins[i] = m.x
            values[i] = m.fun

        i1 = np.where(values == values.min())
        values = np.delete(values, i1)
        i2 = np.where(values == values.min())

        y1 = mins[i1]
        y2 = mins[i2]
        x1 = c[0, i1]
        x2 = c[0, i2]

        coords_min1 = np.array([x1, y1])
        coords_min2 = np.array([x2, y2])

    return coords_min1, coords_min2

coords_min1, coords_min2 = evaluate_minimum(dimension='y')

i = Inflow(t_path, snapshot)
i.read_table()

plt.plot(i.x, i.y, ',')
plt.plot(coords_min1[0], coords_min1[1], 'X')
plt.plot(coords_min2[0], coords_min2[1], 'X')
plt.show()
