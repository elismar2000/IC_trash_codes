from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from gas_to_the_center import r_cm
from particle_positions import position_all_particles, position


def ke(table_path, snapshot, galaxy):
    t = Table.read(table_path, snapshot)
    galaxy_mask = t['GAL'] == galaxy
    v_x = t['VX'][galaxy_mask]
    v_y = t['VY'][galaxy_mask]
    v_z = t['VZ'][galaxy_mask]
    mass = t['MASS'][galaxy_mask]
    ke = mass*(v_x**2 + v_y**2 + v_z**2)/2
    return ke

# def cm(table_path, snapshot, galaxy, energy_threshold):
#     t = Table.read(table_path, snapshot)
#     galaxy_mask = t['GAL'] == galaxy
#     v_x = t['VX'][galaxy_mask]
#     v_y = t['VY'][galaxy_mask]
#     v_z = t['VZ'][galaxy_mask]
#     mass = t['MASS'][galaxy_mask]
#     kinectic_energy = ke(table_path, snapshot, galaxy)
#     ke_mask = kinectic_energy < energy_threshold
#     x = t['X'][galaxy_mask][ke_mask]
#     y = t['Y'][galaxy_mask][ke_mask]
#     z = t['Z'][galaxy_mask][ke_mask]
#     m = t['MASS'][galaxy_mask][ke_mask]
#     cm = r_cm(x, y, z, m)
#     return cm

if __name__ == '__main__':

    snapshot = 50
    galaxy = 2
    energy_threshold = 200

    table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"

    positions = position_all_particles(table_path, snapshot=snapshot, galaxy=galaxy)
    gas_from_spiral = position(table_path, snapshot=snapshot, galaxy=galaxy, p_type=0)

    ke = ke(table_path, snapshot=snapshot, galaxy=galaxy)

    # cm = cm(table_path, snapshot=snapshot, galaxy=galaxy, energy_threshold=energy_threshold)

    plt.plot(positions[0,:], positions[1,:], weights=ke)
    # plt.plot(positions[0,:], positions[1,:],'.',label='Gas from spiral galaxy')
    # plt.plot(gas_from_spiral[0,:], gas_from_spiral[1,:],'.',label='Gas from spiral galaxy')
    # plt.plot(cm[1], cm[2],'D',label='spiral galaxy cm')
    plt.title('snapshot ' + str(snapshot))
    plt.xlabel('Y (KPc)')
    plt.ylabel('Z (Kpc)')
    plt.legend()
    plt.show()
