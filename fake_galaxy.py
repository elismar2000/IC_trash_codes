import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column

mean = [1, 2]
cov = [[0.1, 0], [0, 0.1]]

x, y = np.random.multivariate_normal(mean, cov, 20000).T

x1 = np.ones(1000)*0.002
x2 = np.ones(2000)*0.001
x3 = np.ones(17000)*0.0005
mass = np.concatenate((x1, x2, x3))

output = Table([
    Column(data=x, name='X'),
    Column(data=y, name='Y'),
    #Column(data=z, name='Z'),
    Column(data=mass, name='MASS')
])
output.write('fake_galaxy_3d.fits', overwrite=True)

plt.hist2d(x, y, bins=50, weights=mass)
plt.colorbar()
plt.show()
