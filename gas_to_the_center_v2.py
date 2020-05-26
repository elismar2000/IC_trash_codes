from astropy.table import Table, Column
import numpy as np
import matplotlib.pyplot as plt
from cm_distance import r_cm_bp

number_of_snapshots = 71
# number_of_particles = len(snapshots['snapshot1'])
snapshots_array = np.arange(0,number_of_snapshots,1)*50

def r_cm(x,y,z,m):
    mass = np.sum(m)
    x = np.sum(x*m)/mass
    y = np.sum(y*m)/mass
    z = np.sum(z*m)/mass
    return np.array([x, y, z])


def distance(table_name, snapshot, galaxy):
    t = Table.read(table_name, snapshot)
    galaxy_mask = t['GAL'] == galaxy
    x = np.array([t[i][galaxy_mask] for i in ['X', 'Y', 'Z']])
    cm = r_cm(x[0], x[1], x[2], t['MASS'][galaxy_mask])
    d = np.sqrt(np.sum(np.square(np.transpose(x) - cm),axis=1))
    return d


def mass_inside(table_name, snapshot, galaxy, r_min, p_type):
    t = Table.read(table_name, snapshot)
    galaxy_mask = t['GAL'] == galaxy
    d = distance(table_name, snapshot, galaxy)
    mask = (t['P_TYPE'][galaxy_mask] == p_type) & (d <= r_min)
    mass = np.sum(t['MASS'][galaxy_mask][mask])
    return mass


def gas_mass(table_name, snapshot, galaxy, r_min):
    t = Table.read(table_name, snapshot)
    galaxy_mask = t['GAL'] == galaxy
    d = distance(table_name, snapshot, galaxy)
    mask = (t['P_TYPE'][galaxy_mask] == 0) & (d <= r_min)
    gas_mass = np.sum(t['M_GAS'][galaxy_mask][mask])
    return gas_mass

def formed_star_mass(table_name, snapshot, galaxy, r_min):
    t = Table.read(table_name, snapshot)
    galaxy_mask = t['GAL'] == galaxy
    d = distance(table_name, snapshot, galaxy)
    mask = (t['P_TYPE'][galaxy_mask] == 0) & (d <= r_min)
    m_gas = t['M_GAS'][galaxy_mask][mask]
    mass = t['MASS'][galaxy_mask][mask]
    star_mass = np.sum(mass - m_gas)
    return star_mass


if __name__ == '__main__':
    tn = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"

    if True:
        f = 2.25e9
        dark_matter_galaxy2 = np.array([mass_inside(tn, snapshot=i, galaxy=2, r_min=1, p_type=2) for i in range(1, number_of_snapshots+1)])*f
        dark_matter_galaxy1 = np.array([mass_inside(tn, snapshot=i, galaxy=1, r_min=1, p_type=2) for i in range(1, number_of_snapshots+1)])*f
        gas_galaxy2 = np.array([gas_mass(tn, snapshot=i, galaxy=2, r_min=1) for i in range(1, number_of_snapshots+1)])*f
        original_star_mass_galaxy2 = np.array([mass_inside(tn, snapshot=i, galaxy=2, r_min=1, p_type=1) for i in range(1, number_of_snapshots+1)])
        formed_star_mass_galaxy2 =  np.array([star_mass(tn, snapshot=i, galaxy=2, r_min=1) for i in range(1, number_of_snapshots+1)])
        star_galaxy2 = (original_star_mass_galaxy2 + formed_star_mass_galaxy2)*f
        star_galaxy1 = np.array([mass_inside(tn, snapshot=i, galaxy=1, r_min=1, p_type=1) for i in range(1, number_of_snapshots+1)])

        output = Table([
            Column(data=snapshots_array, name='time'),
            Column(data=dark_matter_galaxy2, name='dark_matter2'),
            Column(data=dark_matter_galaxy1, name='dark_matter1'),
            Column(data=gas_galaxy2, name='gas2'),
            Column(data=star_galaxy2, name='star2'),
            Column(data=star_galaxy1, name='star1'),
        ])
        output.write('output.fits')
    else:
        table = Table.read('output.fits')
        plt.plot(table['time'],table['dark_matter2'],label="Dark matter (galaxy 2)")
        plt.plot(table['time'],table['gas2'],label="Gas (galaxy 2)")
        plt.plot(table['time'],table['star2'],label="Star (galaxy 2)")
        plt.xlabel("t (Myr)")
        plt.ylabel(r"mass ($M_\odot$)")
        plt.legend()
        plt.show()
#
# plt.plot(snapshots_array,dm_mass_array,label='Dark matter')
# plt.plot(snapshots_array,star_mass_array,label='Star')
# plt.plot(snapshots_array,hybrid_mass_array,label='Hybrid')
# plt.legend()
# plt.ylabel('Mass (Msolar)')
# plt.xlabel('Snapshot')
# plt.title('Mass of each type of particle inside a radius of 1 KPc from the system\'s center of mass')
# plt.show()
