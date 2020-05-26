from astropy.table import Table, Column
import numpy as np
import matplotlib.pyplot as plt

def r_cm(x,y,z,m):
    mass = np.sum(m)
    x = np.sum(x*m)/mass
    y = np.sum(y*m)/mass
    z = np.sum(z*m)/mass
    return np.array([x, y, z])


def r_cm_bp(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    particle_mask = table['P_TYPE'][galaxy_mask] != 2
    m = table['MASS'][galaxy_mask][particle_mask]
    x = table['X'][galaxy_mask][particle_mask]
    y = table['Y'][galaxy_mask][particle_mask]
    z = table['Z'][galaxy_mask][particle_mask]
    cm = r_cm(x,y,z,m)
    return cm


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


def distance(table_path, snapshot, galaxy_for_particle, galaxy_for_cm):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy_for_particle
    x = np.array([table[i][galaxy_mask] for i in ['X', 'Y', 'Z']])
    cm = real_cms(table_path, snapshot, galaxy_for_cm)
    d = np.sqrt(np.sum(np.square(np.transpose(x) - cm),axis=1))
    return d


def dark_matter_mass(table_path, snapshot, galaxy_for_particle, galaxy_for_cm):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy_for_particle
    d = distance(table_path, snapshot, galaxy_for_particle, galaxy_for_cm)
    r_min = 1
    mask = (table['P_TYPE'][galaxy_mask] == 2) & (d <= r_min)
    mass = np.sum(table['MASS'][galaxy_mask][mask])
    return mass


def gas_mass(table_path, snapshot, galaxy_for_particle, galaxy_for_cm):
    t = Table.read(table_path, snapshot)
    galaxy_mask = t['GAL'] == galaxy_for_particle
    r_min = 1
    d = distance(table_path, snapshot, galaxy_for_particle, galaxy_for_cm)
    mask = (t['P_TYPE'][galaxy_mask] == 0) & (d <= r_min)
    gas_mass = np.sum(t['M_GAS'][galaxy_mask][mask])
    return gas_mass


def formed_star_mass(table_path, snapshot, galaxy_for_particle, galaxy_for_cm):
    t = Table.read(table_path, snapshot)
    galaxy_mask = t['GAL'] == galaxy_for_particle
    d = distance(table_path, snapshot, galaxy_for_particle, galaxy_for_cm)
    r_min = 1
    gas_mask = (t['P_TYPE'][galaxy_mask] == 0) & (d <= r_min)
    m_gas = t['M_GAS'][galaxy_mask][gas_mask]
    mass = t['MASS'][galaxy_mask][gas_mask]
    formed_star_mass = np.sum(mass - m_gas)
    return formed_star_mass


def star_mass(table_path, snapshot, galaxy_for_particle, galaxy_for_cm):
    t = Table.read(table_path, snapshot)
    galaxy_mask = t['GAL'] == galaxy_for_particle
    d = distance(table_path, snapshot, galaxy_for_particle, galaxy_for_cm)
    r_min = 1
    star_mask = (t['P_TYPE'][galaxy_mask] == 1) & (d <= r_min)
    star_mass = t['MASS'][galaxy_mask][star_mask]
    formed_star_mass_ = formed_star_mass(table_path, snapshot, galaxy_for_particle, galaxy_for_cm)
    total_star_mass = np.sum(star_mass + formed_star_mass_)
    return total_star_mass


if __name__ == '__main__':

    table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
    number_of_snapshots = 71
    snapshots_array = np.arange(0,number_of_snapshots,1)*50

    un_massa = 2.25e9

    # dm_galaxy1_fromcm1 = np.array([dark_matter_mass(table_path, i, 1, 1) for i in range(1, number_of_snapshots+1)])*un_massa
    # dm_galaxy1_fromcm2 = np.array([dark_matter_mass(table_path, i, 1, 2) for i in range(1, number_of_snapshots+1)])*un_massa
    # dm_galaxy1 = dm_galaxy1_fromcm1 + dm_galaxy1_fromcm2

    gas_galaxy2_fromcm1 = np.array([gas_mass(table_path, i, 2, 1) for i in range(1, number_of_snapshots+1)])*un_massa
    gas_galaxy2_fromcm2 = np.array([gas_mass(table_path, i, 2, 2) for i in range(1, number_of_snapshots+1)])*un_massa
    gas_galaxy2 = gas_galaxy2_fromcm1 + gas_galaxy2_fromcm2

    # star_galaxy2_fromcm1 = np.array([star_mass(table_path, i, 2, 1) for i in range(1, number_of_snapshots+1)])*un_massa
    # star_galaxy2_fromcm2 = np.array([star_mass(table_path, i, 2, 2) for i in range(1, number_of_snapshots+1)])*un_massa
    # star_galaxy2 = star_galaxy2_fromcm1 + star_galaxy2_fromcm2

    # star_galaxy1_fromcm1 = np.array([star_mass(table_path, i, 1, 1) for i in range(1, number_of_snapshots+1)])*un_massa
    # star_galaxy1_fromcm2 = np.array([star_mass(table_path, i, 1, 2) for i in range(1, number_of_snapshots+1)])*un_massa
    # star_galaxy1 = star_galaxy1_fromcm1 + star_galaxy1_fromcm2

    plt.plot(snapshots_array, gas_galaxy2_fromcm1, label='cm1')
    plt.plot(snapshots_array, gas_galaxy2_fromcm2, label='cm2')
    plt.plot(snapshots_array, gas_galaxy2, label='sum')
    plt.xlabel('time (Myr)')
    plt.ylabel(r'mass ($M_\odot$)')
    plt.title('Mass of gas particles from spiral galaxy for r_min = 1 Kpc')
    plt.legend()
    plt.show()
