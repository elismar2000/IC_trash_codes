import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import Animation


def limits(table_path, snapshot, bins, delta, galaxy):
    t = Table.read(table_path, snapshot)
    g_mask = t['GAL'] == galaxy
    h1, xedges, yedges = np.histogram2d(t['X'][g_mask], t['Y'][g_mask], bins=[bins, bins], weights=t['MASS'][g_mask])
    h2, aux1, zedges = np.histogram2d(t['X'][g_mask], t['Z'][g_mask], bins=[bins, bins], weights=t['MASS'][g_mask])
    a, b = np.where(h1 == np.amax(h1))
    aux2, c = np.where(h2 == np.amax(h2))
    x_llim = xedges[a[0] - delta]
    x_ulim = xedges[a[0] + delta]
    y_llim = yedges[b[0] - delta]
    y_ulim = yedges[b[0] + delta]
    z_llim = zedges[c[0] - delta]
    z_ulim = zedges[c[0] + delta]
    limits = np.array([x_llim, x_ulim, y_llim, y_ulim, z_llim, z_ulim])
    return limits


def bins_mask(table_path, snapshot, bins, delta):
    t = Table.read(table_path, snapshot)
    limits1 = limits(table_path, snapshot, bins, delta, galaxy=1)
    limits2 = limits(table_path, snapshot, bins, delta, galaxy=2)
    x_llim = np.amin(np.array([limits1[0], limits2[0]]))
    x_ulim = np.amax(np.array([limits1[1], limits2[1]]))
    y_llim = np.amin(np.array([limits1[2], limits2[2]]))
    y_ulim = np.amax(np.array([limits1[3], limits2[3]]))
    z_llim = np.amin(np.array([limits1[4], limits2[4]]))
    z_ulim = np.amax(np.array([limits1[5], limits2[5]]))
    mask = ((x_llim < t['X']) & (t['X'] < x_ulim)) & ((y_llim < t['Y']) & (t['Y'] < y_ulim)) & ((z_llim < t['Z']) & (t['Z'] < z_ulim))
    return mask


def limits_again(tale_path, snapshot, bins, delta1, delta2, galaxy):
    t = Table.read(table_path, snapshot)
    mask = bins_mask(table_path, snapshot, bins, delta1)
    g_mask = t['GAL'] == galaxy
    mask = mask[g_mask]
    h1, xedges, yedges = np.histogram2d(t['X'][g_mask][mask], t['Y'][g_mask][mask], bins=[bins, bins], weights=t['MASS'][g_mask][mask])
    h2, aux1, zedges = np.histogram2d(t['X'][g_mask][mask], t['Z'][g_mask][mask], bins=[bins, bins], weights=t['MASS'][g_mask][mask])
    a, b = np.where(h1 == np.amax(h1))
    aux2, c = np.where(h2 == np.amax(h2))
    x_llim = xedges[a[0] - delta2]
    x_ulim = xedges[a[0] + delta2]
    y_llim = yedges[b[0] - delta2]
    y_ulim = yedges[b[0] + delta2]
    z_llim = zedges[c[0] - delta2]
    z_ulim = zedges[c[0] + delta2]
    limits = np.array([x_llim, x_ulim, y_llim, y_ulim, z_llim, z_ulim])
    return limits


def bins_mask_again(table_path, snapshot, bins, delta1, delta2):
    t = Table.read(table_path, snapshot)
    limits1 = limits_again(table_path, snapshot, bins, delta1, delta2, galaxy=1)
    limits2 = limits_again(table_path, snapshot, bins, delta1, delta2, galaxy=2)
    x_llim = np.amin(np.array([limits1[0], limits2[0]]))
    x_ulim = np.amax(np.array([limits1[1], limits2[1]]))
    y_llim = np.amin(np.array([limits1[2], limits2[2]]))
    y_ulim = np.amax(np.array([limits1[3], limits2[3]]))
    z_llim = np.amin(np.array([limits1[4], limits2[4]]))
    z_ulim = np.amax(np.array([limits1[5], limits2[5]]))
    mask = ((x_llim < t['X']) & (t['X'] < x_ulim)) & ((y_llim < t['Y']) & (t['Y'] < y_ulim)) & ((z_llim < t['Z']) & (t['Z'] < z_ulim))
    return mask


if __name__ == '__main__':
    table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
    snapshot = 50
    bins = 200
    delta1 = 2
    delta2 = 4
    mask = bins_mask_again(table_path, snapshot=snapshot, bins=bins, delta1=delta1, delta2=delta2)

    t = Table.read(table_path, snapshot)
    fig, ax = plt.subplots(1,3, figsize=(10,10), tight_layout=True)
    hist1 = ax[0].hist2d(t['X'][mask], t['Y'][mask], bins=(bins, bins), weights=t['MASS'][mask])
    ax[0].set_aspect('equal')

    hist2 = ax[1].hist2d(t['X'][mask], t['Z'][mask], bins=(bins, bins), weights=t['MASS'][mask])
    ax[1].set_aspect('equal')

    hist3 = ax[2].hist2d(t['Y'][mask], t['Z'][mask], bins=(bins, bins), weights=t['MASS'][mask])
    ax[2].set_aspect('equal')

    fig.tight_layout()
    plt.show()
