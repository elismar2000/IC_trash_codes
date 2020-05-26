from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)

def galstep(r, mass, a):
    return (mass / (2 * np.pi)) * (a / (r * (r + a)**3))

gmodel = Model(gaussian)
galmodel = Model(galstep)

gauss_params = gmodel.make_params(cen=5, amp=200, wid=1)
galstep_params = galmodel.make_params(mass=2.3e+10, a=1.5)

x = np.linspace(0.1, 10, 201)
#y_gauss = gmodel.eval(gauss_params, x=x)
#y_galstep = galmodel.eval(galstep_params, r=x)
y_gauss = gmodel.eval(x=x, cen=6.5, amp=100, wid=2.0)
y_galstep = galmodel.eval(r=x, mass=3.0e+10, a=3.0)

result_gauss = gmodel.fit(y_gauss, gauss_params, x=x)
result_galstep = galmodel.fit(y_galstep, galstep_params, r=x)

print(result_gauss.fit_report())
print(result_galstep.fit_report())

fig, axs = plt.subplots(1, 2)
axs[0].plot(x, y_gauss, 'o')
axs[0].plot(x, result_gauss.init_fit, 'k--', label='init fit')
axs[0].plot(x, result_gauss.best_fit, 'r-', label='best fit')
axs[0].legend()

axs[1].plot(x, y_galstep, 'o')
axs[1].plot(x, result_galstep.init_fit, 'k--', label='init fit')
axs[1].plot(x, result_galstep.best_fit, 'r-', label='best fit')
axs[1].legend()

plt.show()
