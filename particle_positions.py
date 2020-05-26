from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from cm_distance import r_cm_ap, r_cm_bp, r_cm_wsbp, r_cm_wsap
#from potential import MyTable


def position(table_path, snapshot, galaxy, p_type):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    particle_mask = table['P_TYPE'][galaxy_mask] == p_type
    x = table['X'][galaxy_mask][particle_mask]
    z = table['Z'][galaxy_mask][particle_mask]
    y = table['Y'][galaxy_mask][particle_mask]
    return np.array([x, y, z])


def position_both_galaxies(table_path, snapshot, p_type):
    table = Table.read(table_path, snapshot)
    particle_mask = table['P_TYPE'] == p_type
    x = table['X'][particle_mask]
    z = table['Z'][particle_mask]
    y = table['Y'][particle_mask]
    return np.array([x, y, z])


def position_all_particles(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    x = table['X'][galaxy_mask]
    y = table['Y'][galaxy_mask]
    z = table['Z'][galaxy_mask]
    return np.array([x, y, z])


def position_all(table_path, snapshot):
    table = Table.read(table_path, snapshot)
    x = table['X']
    y = table['Y']
    return np.array([x, y])


def r_cm(x,y,z,m):
    mass = np.sum(m)
    x = np.sum(x*m)/mass
    y = np.sum(y*m)/mass
    z = np.sum(z*m)/mass
    return np.array([x, y, z])


def select_particles(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    particle_mask = table['P_TYPE'] != 2
    x = table['X'][particle_mask]
    y = table['Y'][particle_mask]
    z = table['Z'][particle_mask]
    r = np.array([x,y,z])
    radius_for_cm = 33
    cm_bp = r_cm_bp(table_path, snapshot, galaxy)
    mask = np.sqrt(np.sum(np.square(np.transpose(r) - cm_bp),axis=1)) < radius_for_cm
    return mask


def real_cms(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    mask = select_particles(table_path, snapshot, galaxy)
    particle_mask = table['P_TYPE'] != 2
    x = table['X'][particle_mask][mask]
    y = table['Y'][particle_mask][mask]
    z = table['Z'][particle_mask][mask]
    m = table['MASS'][particle_mask][mask]
    cm = r_cm(x,y,z,m)
    return cm

if __name__ == '__main__':

    table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
    snapshot = 71
    galaxy1 = 1
    galaxy2 = 2
    n_bins = 20

#    t1 = MyTable(table_path, snapshot, galaxy1, n_bins = n_bins)
#    t2 = MyTable(table_path, snapshot, galaxy2, n_bins = n_bins)
#    u1, coords_min1 = t1.potential()
#    u2, coords_min2 = t2.potential()

#    r_cm_bp1 = r_cm_bp(table_path, snapshot, galaxy1)
#    r_cm_bp2 = r_cm_bp(table_path, snapshot, galaxy2)

    number_of_snapshots = 71
    snapshots_array = np.arange(0,number_of_snapshots,1)*0.5

    positions1 = position_all_particles(table_path, snapshot=snapshot, galaxy=galaxy1)
    positions2 = position_all_particles(table_path, snapshot=snapshot, galaxy=galaxy2)
    gas_from_spiral = position(table_path, snapshot=snapshot, galaxy=2, p_type=0)
    dm = position_both_galaxies(table_path, snapshot, 2)

#    coords1 = t1.coordinates()
#    lado_do_bin = coords1[0, 1, 1, 1] - coords1[0, 1, 0, 1]

    plt.plot(positions1[0, :], positions1[1, :], ',', label='Particles from eliptical galaxy')
    plt.plot(positions2[0, :], positions2[1, :], ',', label='Particles from spiral galaxy')
    plt.plot(gas_from_spiral[0, :], gas_from_spiral[1, :], ',', label='Gas from spiral galaxy')
#    plt.plot(dm[0, :], dm[1, :], ',', label='dark matter')
    plt.title('snapshot ' + str(snapshot))
    plt.xlabel('X (KPc)')
    plt.ylabel('Y (Kpc)')
    plt.legend()
    plt.show()

    # hexbin(positions1[0,:], positions1[1,:], extent=[-100, 100, -100, 100])
    # plt.plot(positions1[0,:], positions[1,:], '.', label='particles from eliptical galaxy')
    # plt.plot(gas_from_spiral[])
