from scipy.optimize import curve_fit, minimize
import numpy as np
import matplotlib.pyplot as plt

c = np.linspace(1, 10, 1000)
pot = -1 / c**2

#Polinônio de Legrendre P7(x)
def func(x, a, b, c, d):
    return a*x**7 - b*x**5 + c*x**3 - d*x

#Coeficientes ajustados
popt, popv = curve_fit(func, c, pot)

#Reconstruindo a curva pra plotar
p = func(c, popt[0], popt[1], popt[2], popt[3])

plt.plot(c, pot, '--', label='binned 1d potential')
plt.plot(c, p, label='best fit')

#Função que melhor ajusta o potencial
def best_fit(x):
    return func(x, popt[0], popt[1], popt[2], popt[3])

#Minimizando a curva
min = minimize(best_fit, [0], bounds=[[c.min(), c.max()]])

plt.plot(min.x, min.fun, 'X', label='pot min')

#Criando uma curva com muito mais pontos pra mostrar que o mínimo tá no lugar certo
x = np.arange(c.min(), c.max(), 0.000001)
f = best_fit(x)

plt.plot(x, f, label='best fit with much more points')

plt.legend()
plt.show()
