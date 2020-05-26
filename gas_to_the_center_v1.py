from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

table = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
number_of_snapshots = 71

snapshots = {"snapshot" + str(i): Table.read(table,i) for i in range(1,number_of_snapshots+1)}

number_of_particles = len(snapshots['snapshot1'])
snapshots_array = np.arange(1,number_of_snapshots+1,1)

x = np.vstack([snapshots['snapshot'+str(j)]['X']] for j in range(1,number_of_snapshots+1))
y = np.vstack([snapshots['snapshot'+str(j)]['Y']] for j in range(1,number_of_snapshots+1))
z = np.vstack([snapshots['snapshot'+str(j)]['Z']] for j in range(1,number_of_snapshots+1))

r = np.array([x,y,z])

vx = np.vstack([snapshots['snapshot'+str(j)]['VX']] for j in range(1,number_of_snapshots+1))
vy = np.vstack([snapshots['snapshot'+str(j)]['VY']] for j in range(1,number_of_snapshots+1))
vz = np.vstack([snapshots['snapshot'+str(j)]['VZ']] for j in range(1,number_of_snapshots+1))

v = np.array([vx,vy,vz])

mass = np.vstack([snapshots['snapshot'+str(j)]['MASS']] for j in range(1,number_of_snapshots+1))

x_cm = np.sum(x*mass,axis=1)/np.sum(mass,axis=1)
y_cm = np.sum(y*mass,axis=1)/np.sum(mass,axis=1)
z_cm = np.sum(z*mass,axis=1)/np.sum(mass,axis=1)

r_cm = np.array([x_cm,y_cm,z_cm])

r_radial = np.zeros_like(r)
for j in range(0,number_of_particles):
    r_radial[:,:,j] = r[:,:,j] - r_cm

r_normalized = r_radial/np.sqrt(np.sum(r_radial*r_radial,axis=0))
gas_to_the_center = np.sum(v*r_normalized,axis=0)
the_last = np.sum(gas_to_the_center,axis=1)

plt.plot(snapshots_array,the_last,label='')
plt.ylabel('Flux of particles to the center of the galaxies')
plt.xlabel('Snapshot')
plt.show()

# x = np.vstack([snapshots['snapshot'+str(j)]['X']] for j in range(1,number_of_snapshots+1))
# y = np.vstack([snapshots['snapshot'+str(j)]['Y']] for j in range(1,number_of_snapshots+1))
# z = np.vstack([snapshots['snapshot'+str(j)]['Z']] for j in range(1,number_of_snapshots+1))
#
# p_type = np.vstack([snapshots['snapshot'+str(j)]['P_TYPE']] for j in range(1,number_of_snapshots+1))
#
# r = np.array([x,y,z])
#
# mass = np.vstack([snapshots['snapshot'+str(j)]['MASS']] for j in range(1,number_of_snapshots+1))
#
# r_cm = np.array([x_cm,y_cm,z_cm])
#
# r_radial = np.zeros_like(r)
# for j in range(0,number_of_particles):
#     r_radial[:,:,j] = r[:,:,j] - r_cm
#
# r_length = np.sqrt(np.sum(r_radial*r_radial,axis=0))
#
# dm_mass = {'mass' + str(i): mass[i,(r_length[i,:] <= 1) & (p_type[i,:] == 2)] for i in range(0,number_of_snapshots)}
# star_mass = {'mass' + str(i): mass[i,(r_length[i,:] <= 1) & (p_type[i,:] == 1)] for i in range(0,number_of_snapshots)}
# hybrid_mass = {'mass' + str(i): mass[i,(r_length[i,:] <= 1) & (p_type[i,:] == 0)] for i in range(0,number_of_snapshots)}
#
# dm_mass_array = np.array([np.sum(dm_mass['mass'+str(i)]) for i in range(0,number_of_snapshots)])
# star_mass_array = np.array([np.sum(star_mass['mass'+str(i)]) for i in range(0,number_of_snapshots)])
# hybrid_mass_array = np.array([np.sum(hybrid_mass['mass'+str(i)]) for i in range(0,number_of_snapshots)])
