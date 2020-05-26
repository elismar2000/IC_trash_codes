import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.table import Table

class Potential:
    def __init__(self):
        self.coords = np.array([])
        self.mass = np.array([])
        self.mask = np.array([])
        self.d_bin = 0.0

    def read_table(self, t_path):
        t = Table.read(t_path)
        self.coords = t['C']
        self.mass = t['MASS']

    def evaluate_potential(self, n_bins: int = 10, smooth: float = 2.0):
        h, bins = np.histogram(self.coords, bins=n_bins, weights=self.mass)
        n, bins2 = np.histogram(self.coords, bins=n_bins)
        d_bin = (bins[-1] - bins[0]) / n_bins
        self.d_bin = d_bin
        centers = bins[:-1] + (d_bin / 2.0)

        potential = np.zeros_like(h)
        for i, u in enumerate(potential):
            distance = np.abs(centers - centers[i])
            distance[i] = d_bin / 4.0
            potential += -(h / distance)

        if smooth != 0.0:
            potential = convolve(potential, Gaussian1DKernel(smooth))

        return centers, potential, n

    def select(self, x0: float, width: float = 1.5):
        w = self.d_bin * width
        mask = np.abs(self.coords - x0) < w
        self.coords = self.coords[mask]
        self.mass = self.mass[mask]

    def pot1d_2(self, t_path, n_bins: int = 10, smooth: float = 1.5, width: float = 2.0):
        self.read_table(t_path)
        c, p, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)
        d = np.diff(p)
        zero = np.array([(d[i] > 0) & (d[i-1] < 0) for i in range(len(d))])
        [i] = np.where(zero == True)
        x_min1 = c[i[0]]
        x_min2 = c[i[1]]

        x1 = self.converge(x_min1, t_path, n_bins, smooth, width)
        x2 = self.converge(x_min2, t_path, n_bins, smooth, width)

        return x1, x2

    def converge(self, x, t_path, n_bins: int = 10, smooth: float = 1.5, width: float = 1.5):
        n_min = 10.0
        self.read_table(t_path=t_path)
        self.evaluate_potential(n_bins=n_bins, smooth=smooth)
        self.select(x, width=width)
        while n_min > 2.0:
            c, p, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)
            d = np.diff(p)
            zero = np.array([(d[i] > 0) & (d[i-1] < 0) for i in range(len(d))])
            [i] = np.where(zero == True)
            x_min = c[i[0]]

            self.select(x_min, width=width)

            n_min = n.min()
            print(x_min, n_min)

        return x_min


if __name__ == '__main__':

    p = Potential()

    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/fake_galaxy.fits'

    x1, x2 = p.pot1d_2(t_path=t_path, n_bins=10, smooth=2.0, width=2.0)
