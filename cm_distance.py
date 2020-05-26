from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
number_of_snapshots = 71

def r_cm_ap(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    m = table['MASS'][galaxy_mask]
    x = table['X'][galaxy_mask]
    y = table['Y'][galaxy_mask]
    z = table['Z'][galaxy_mask]
    mass = np.sum(m)
    x_cm = np.sum(x*m)/mass
    y_cm = np.sum(y*m)/mass
    z_cm = np.sum(z*m)/mass
    return np.array([x_cm, y_cm, z_cm])


def r_cm_bp(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    particle_mask = table['P_TYPE'][galaxy_mask] != 2
    m = table['MASS'][galaxy_mask][particle_mask]
    x = table['X'][galaxy_mask][particle_mask]
    y = table['Y'][galaxy_mask][particle_mask]
    z = table['Z'][galaxy_mask][particle_mask]
    mass = np.sum(m)
    x_cm = np.sum(x*m)/mass
    y_cm = np.sum(y*m)/mass
    z_cm = np.sum(z*m)/mass
    return np.array([x_cm, y_cm, z_cm])


def r_cm_wsap(table_path, snapshot):
    table = Table.read(table_path, snapshot)
    m = table['MASS']
    x = table['X']
    y = table['Y']
    z = table['Z']
    mass = np.sum(m)
    x_cm = np.sum(x*m)/mass
    y_cm = np.sum(y*m)/mass
    z_cm = np.sum(z*m)/mass
    return np.array([x_cm, y_cm, z_cm])


def r_cm_wsbp(table_path, snapshot):
    table = Table.read(table_path, snapshot)
    particle_mask = table['P_TYPE'] != 2
    m = table['MASS'][particle_mask]
    x = table['X'][particle_mask]
    y = table['Y'][particle_mask]
    z = table['Z'][particle_mask]
    mass = np.sum(m)
    x_cm = np.sum(x*m)/mass
    y_cm = np.sum(y*m)/mass
    z_cm = np.sum(z*m)/mass
    return np.array([x_cm, y_cm, z_cm])


def r_cm_distance(cm1,cm2):
    r_cm_distance_vector = cm1 - cm2
    r_cm_distance = np.sqrt(np.sum(np.square(r_cm_distance_vector)))
    return r_cm_distance

if __name__ == '__main__':
    cm1_ap = np.array([r_cm_ap(table_path, i, 1) for i in range(1, number_of_snapshots+1)])
    cm2_ap = np.array([r_cm_ap(table_path, i, 2) for i in range(1, number_of_snapshots+1)])
    cm1_bp = np.array([r_cm_bp(table_path, i, 1) for i in range(1, number_of_snapshots+1)])
    cm2_bp = np.array([r_cm_bp(table_path, i, 2) for i in range(1, number_of_snapshots+1)])
    cm_wsap = np.array([r_cm_wsap(table_path, i) for i in range(1, number_of_snapshots+1)])
    cm_wsbp = np.array([r_cm_wsbp(table_path, i) for i in range(1, number_of_snapshots+1)])
    r_cm_distance_barionic = np.array([r_cm_distance(cm1_bp[i], cm2_bp[i]) for i in range(0, len(cm1_bp))])
    r_cm_distance_all = np.array([r_cm_distance(cm1_ap[i], cm2_ap[i]) for i in range(0, len(cm1_ap))])

    snapshots_array = np.arange(0,number_of_snapshots,1)*50
    # plt.plot(cm1_ap[:,0], cm1_ap[:,1],label='eliptical galaxy - all particles')
    # plt.plot(cm2_ap[:,0], cm2_ap[:,1],label='spiral galaxy - all particles')
    # plt.plot(cm1_bp[:,0], cm1_bp[:,1],label='eliptical galaxy - barionic particles')
    # plt.plot(cm2_bp[:,0], cm2_bp[:,1],label='spiral galaxy - barionic particles')
    # plt.plot(cm_wsap[:,0], cm_wsap[:,1],label='whole system - all particles')
    # plt.plot(cm_wsbp[:,0], cm_wsbp[:,1],label='whole system - barionic particles')

    plt.plot(snapshots_array, r_cm_distance_barionic, label='barionic')
    plt.plot(snapshots_array, r_cm_distance_all, label='all')
    plt.xlabel('Tempo decorrido do início a colisão (Myr)')
    plt.ylabel('Distância entre os cms em KPc')
    plt.legend()
    plt.show()
