from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from cm_distance import r_cm_bp, r_cm_wsbp
from particle_positions import position_all_particles, position, position_all
from gas_to_the_center_v2 import r_cm

def select_particles_given_galaxy(table_path, snapshot, galaxy, radius, cm):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    x = table['X'][galaxy_mask]
    y = table['Y'][galaxy_mask]
    z = table['Z'][galaxy_mask]
    r = np.array([x,y,z])
    r_cm = cm[snapshot-1,:]
    distance_mask = np.sqrt(np.sum(np.square(np.transpose(r) - r_cm),axis=1)) < radius
    return distance_mask


def select_all_particles(table_path, snapshot, radius, cm):
    table = Table.read(table_path, snapshot)
    x = table['X']
    y = table['Y']
    z = table['Z']
    r = np.array([x,y,z])
    r_cm = cm[snapshot-1,:]
    distance_mask = np.sqrt(np.sum(np.square(np.transpose(r) - r_cm),axis=1)) < radius
    return distance_mask


def real_cm(table_path, snapshot, mask):
    table = Table.read(table_path, snapshot)
    x = table['X'][mask[snapshot-1,:]]
    y = table['Y'][mask[snapshot-1,:]]
    z = table['Z'][mask[snapshot-1,:]]
    m = table['MASS'][mask[snapshot-1,:]]
    cm = r_cm(x,y,z,m)
    return cm


if __name__ == '__main__':
    table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
    number_of_snapshots = 71
    snapshot = 60
    radius = 15

    cm1_bp = np.array([r_cm_bp(table_path, i, galaxy=1) for i in range(1, number_of_snapshots+1)])
    cm2_bp = np.array([r_cm_bp(table_path, i, galaxy=2) for i in range(1, number_of_snapshots+1)])

    mask1 = np.array([select_all_particles(table_path, i, radius=radius, cm=cm1_bp) for i in range(1, number_of_snapshots+1)])
    mask2 = np.array([select_all_particles(table_path, i, radius=radius, cm=cm2_bp) for i in range(1, number_of_snapshots+1)])

    real_cm1 = np.array([real_cm(table_path, i, mask1) for i in range(1, number_of_snapshots+1)])
    real_cm2 = np.array([real_cm(table_path, i, mask2) for i in range(1, number_of_snapshots+1)])

#    positions = position_all(table_path, snapshot=snapshot)
    eliptical_positions = position_all_particles(table_path, snapshot=snapshot, galaxy=1)
    gas_from_spiral = position(table_path, snapshot=snapshot, galaxy=2, p_type=0)

    plt.plot(eliptical_positions[0,:], eliptical_positions[1,:],'.', label='particles from both galaxies')
    plt.plot(gas_from_spiral[0,:], gas_from_spiral[1,:],'.',label='Gas from spiral galaxy')
    # plt.plot(x1,y1,'.', label='particles inside radius of ' + str(radius) + ' kpc from cm1_bp')
    # plt.plot(x2,y2,'.', label='particles inside radius of ' + str(radius) + ' kpc from cm2_bp')
    plt.plot(cm1_bp[snapshot-1,0],cm1_bp[snapshot-1,1],'D', label='cm1_bp')
    plt.plot(cm2_bp[snapshot-1,0],cm2_bp[snapshot-1,1],'D', label='cm2_bp')
    plt.plot(real_cm1[snapshot-1,0],real_cm1[snapshot-1,1],'D', label='real_cm1')
    plt.plot(real_cm2[snapshot-1,0],real_cm2[snapshot-1,1],'D', label='real_cm2')
    plt.title('snapshot '+str(snapshot))
    plt.legend()
    plt.xlabel('X (KPc)')
    plt.ylabel('Y (KPc)')
    plt.show()

    # mask1 = np.array([select_particles_given_galaxy(table_path, i, galaxy=1, radius=radius, cm=cm1_bp) for i in range(1, number_of_snapshots+1)])
    # mask2 = np.array([select_particles_given_galaxy(table_path, i, galaxy=1, radius=radius, cm=cm2_bp) for i in range(1, number_of_snapshots+1)])

#    eliptical_positions = position_all(table_path, snapshot=snapshot, galaxy=1)
    #
    # x1 = eliptical_positions[0,mask1[snapshot-1,:]]
    # y1 = eliptical_positions[1,mask1[snapshot-1,:]]
    # x2 = eliptical_positions[0,mask2[snapshot-1,:]]
    # y2 = eliptical_positions[1,mask2[snapshot-1,:]]
