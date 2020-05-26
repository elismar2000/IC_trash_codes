import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.table import Table
from inflow import Inflow
from relative_minima import find_minima


class Potential:
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])
        self.mass = np.array([])
        self.mask = np.array([])
        self.dx = 0.0
        self.dy = 0.0

    def read_table(self, t_path, snapshot: int = None):
        if snapshot == None:
            t = Table.read(t_path)
        else:
            t = Table.read(t_path, snapshot)
        self.x = t['X']
        self.y = t['Y']
        self.mass = t['MASS']

    def evaluate_potential(self, n_bins: int = 20, smooth: float = 1.3):
        masses, edges = np.histogramdd(np.column_stack([self.x.data, self.y.data]), bins=n_bins, weights=self.mass)
        edges = np.asarray(edges)

        n, edges2 = np.histogramdd(np.column_stack([self.x.data, self.y.data]), bins=n_bins)

        dx = (edges[0][-1] - edges[0][0]) / n_bins
        dy = (edges[1][-1] - edges[1][0]) / n_bins
        self.dx = dx
        self.dy = dy

        center = np.array([edges[0][:-1] + dx / 2, edges[1][:-1] + dy / 2])
        coords = np.meshgrid(*center)

        potential = np.zeros_like(masses)
        for i, u in np.ndenumerate(potential):
            distance = np.sqrt(np.sum([np.square(coords[j] - coords[j][i]) for j in range(2)], 0))
            distance[i] = dx / 4.0
            potential += -(masses / distance)

        if smooth != 0.0:
            potential = convolve(potential, Gaussian2DKernel(smooth))

        return potential, center, n

    def select(self, coords, width: float = 0.75):
        wx = self.dx * width
        wy = self.dy * width
        mask = (np.abs(self.x - coords[0]) < wx) & (np.abs(self.y - coords[1]) < wy)
        self.x = self.x[mask]
        self.y = self.y[mask]
        self.mass = self.mass[mask]

    def pot2d(self, t_path, snapshot: int = None, n_bins: int = 20, smooth: float = 1.3, width: float = 1.5):
        self.read_table(t_path, snapshot)
        n_min = 200.0
        while n_min > 100.0:
            p, c, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)

            index_min = np.asarray(find_minima(p))
            mins = p[index_min[:, 0], index_min[:, 1]]
            print('mins: ', mins)
            i = np.where(mins < mins.min() / 2)

            j = []
            for _ in range(np.size(i)):
                j.append(np.where(p == mins[i[0][_]]))

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(p, origin='lower')
            for _ in range(np.size(j, axis=0)):
                ax.scatter(j[_][1], j[_][0])
            plt.show()

            if np.size(j, axis=0) == 1:
                coords_min = np.array([float(c[i][j[0][i]]) for i in range(2)])
                self.select(coords_min, width=width)
                n_min = n[j[0]][0]
                print('tinha só um mínimo')

            if np.size(j, axis=0) == 2:
                dx, dy = self.dx, self.dy
                coords_min1 = np.array([float(c[i][j[0][i]]) for i in range(2)])
                coords_min1 = self.converge(coords_min1, t_path, snapshot, n_bins=n_bins, smooth=smooth, width=width)

                self.dx, self.dy = dx, dy
                coords_min2 = np.array([float(c[i][j[1][i]]) for i in range(2)])
                coords_min2 = self.converge(coords_min2, t_path, snapshot, n_bins=n_bins, smooth=smooth, width=width)

                coords_min = np.array([coords_min1, coords_min2])
                print('tinha só dois mínimos')

                break

            if np.size(j, axis=0) > 2:
                i = np.where(p == p.min())
                assert all(len(_) == 1 for _ in i), "Deu mais de um minimo."
                coords_min = np.array([float(c[_][i[_]]) for _ in range(2)])
                self.select(coords_min, width=width)
                n_min = float(n[i])
                print('Tinha mais de dois mínimos')

        return coords_min

    def converge(self, coords_min, t_path, snapshot, n_bins: int = 20, smooth: float = 1.3, width: float = 1.5):
        self.read_table(t_path, snapshot)
        self.select(coords_min, width=width)
        n_min = 200.0
        while n_min > 100.0:
            p, c, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)
            index_min = np.where(p == p.min())
            coords_min = np.array([float(c[i][index_min[i]]) for i in range(2)])
            self.select(coords_min, width=width)

            n_min = n[index_min]

        return coords_min

    def teste(self, t_path, snapshot, n_bins: int = 20, smooth: float = 1.3, width: float = 1.5, zoom: int = 1):
        self.read_table(t_path, snapshot)
        for i in range(zoom):
            p, c, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)
            index_min = np.where(p == p.min())
            coords_min = np.array([float(c[i][index_min[i]]) for i in range(2)])
            self.select(coords_min, width=width)

            n_min = n.min()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(p, origin='lower')
        plt.show()

        return p, c


if __name__ == '__main__':

    p = Potential()

    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit1'
    snapshot = 9

    coords_min = p.pot2d(t_path, snapshot, n_bins=10, smooth=1.3, width=3.0)

    i = Inflow(t_path, snapshot)
    i.read_table()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(i.x, i.y, ',')
    if coords_min.shape == (2,):
        ax.plot(coords_min[0], coords_min[1], 'ko')
    if coords_min.shape == (2, 2):
        ax.plot(coords_min[0, 0], coords_min[0, 1], 'ko')
        ax.plot(coords_min[1, 0], coords_min[1, 1], 'ko')
    plt.show()

    #p.read_table(t_path, snapshot)
    #pot, c, n = p.evaluate_potential(n_bins=20, smooth=1.3)
