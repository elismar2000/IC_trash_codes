from potential_1d import Potential
from astropy.table import Table
import numpy as np
from astropy.modeling import models, fitting, polynomial
import matplotlib.pyplot as plt
from scipy.optimize import minimize

#Calculando o potencial pra uma galáxia gaussiana fake
p = Potential()

coords = np.random.normal(0, 3, 1000)
mass = np.random.choice([1., 2.], 1000)

p.coords = coords
p.mass = mass
c, pot, n = p.evaluate_potential(n_bins=10)

#Ajustando uma gaussiana nessa galáxia com o astropy
g_init = models.Gaussian1D(amplitude=-1., mean=0, stddev=1.)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, c, pot)

#Ajustando um polinômio
p_init = polynomial.Polynomial1D(degree=6)
p = fit_g(p_init, c, pot)

#Minimizando
m = minimize(g, 0)

#plotando
plt.figure(figsize=(8,5))
plt.plot(c, pot, 'ko')
plt.plot(c, g(c), label='best gaussian fit')
plt.plot(c, p(c), label='best sixth degree polynomial fit')
plt.plot(m.x, m.fun, 'X', label='minimum')
plt.legend()
plt.show()
