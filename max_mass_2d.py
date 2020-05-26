import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel
from scipy.ndimage import gaussian_filter
from astropy.table import Table
from particle_positions import position_all_particles, position


class Potential:
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])
        self.mass = np.array([])
        self.mask = np.array([])
        self.dx = 0.0
        self.dy = 0.0

    @staticmethod
    def make_fake_galaxy(table_path, snapshot: int = None):
        if snapshot == None:
            t = Table.read(table_path)
        else:
            t = Table.read(table_path, snapshot)
        x = t['X']
        y = t['Y']
        mass = t['MASS']
        return x, y, mass

    def find_max_mass(self, n_bins: int = 10, smooth: float = 1.0):
        masses, edges = np.histogramdd(np.array([self.x, self.y]).T, bins=n_bins, weights=self.mass)
        edges = np.asarray(edges)

        masses = convolve(masses, Gaussian2DKernel(smooth))

        dx = (edges[0][-1] - edges[0][0]) / n_bins
        dy = (edges[1][-1] - edges[1][0]) / n_bins
        self.dx = dx
        self.dy = dy
        #import pdb; pdb.set_trace()
        index_max = np.where(masses == masses.max())
        x_min = edges[0, index_max[0]]
        y_min = edges[1, index_max[1]]
        x = x_min + dx / 2
        y = y_min + dy / 2

        return masses, edges, x, y

    def select(self, coords, width: float = 0.75):
        wx = self.dx * width
        wy = self.dy * width
        mask = (np.abs(self.x - coords[0]) < wx) & (np.abs(self.y - coords[1]) < wy)
        self.x = self.x[mask]
        self.y = self.y[mask]
        self.mass = self.mass[mask]


if __name__ == '__main__':

    p = Potential()

    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2'
    snapshot = 71
    p.x, p.y, p.mass = p.make_fake_galaxy(t_path, snapshot)

    #m, e, x, y = p.find_max_mass(10)

    positions1 = position_all_particles(t_path, snapshot=snapshot, galaxy=1)
    positions2 = position_all_particles(t_path, snapshot=snapshot, galaxy=2)
    gas_from_spiral = position(t_path, snapshot=snapshot, galaxy=2, p_type=0)

    fig, ax = plt.subplots()
    plt.plot(positions1[0, :], positions1[1, :], ',', label='Particles from eliptical galaxy')
    plt.plot(positions2[0, :], positions2[1, :], ',', label='Particles from spiral galaxy')
    plt.title('snapshot ' + str(snapshot))
    plt.xlabel('X (KPc)')
    plt.ylabel('Y (Kpc)')
    plt.legend()


    for k in range(10):
        m, e, x, y = p.find_max_mass(n_bins=10)
        print('maximum mass: {:.2e}; x_min: {:.2f}, y_min: {:.2f}'.format(m.max(), x[0], y[0]))
        p.select(coords=np.array([x, y]), width=1.0)
        plt.plot(x, y, 'X')
        plt.show()
        input()


# potential = np.zeros_like(masses)
# for i, u in np.ndenumerate(potential):
#     distance = np.sqrt(np.sum([np.square(coords[j] - coords[j][i]) for j in range(2)], 0))
#     distance[i] = dx / 4.0
#     potential += -(masses / distance)
#
# if smooth != 0.0:
#     masses = convolve(masses, Gaussian2DKernel(smooth))
