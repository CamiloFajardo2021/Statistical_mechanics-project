import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

def f(x):
    return 1+x/8

k, df, Ddf = np.loadtxt('data.dat', unpack=True)

lin= stats.linregress(k, df)

plt.plot(k, lin.slope*k+lin.intercept, '--', label='Ajuste')
plt.errorbar(k, df, yerr=Ddf, fmt='ko', capsize=5.0)
plt.plot(k, f(k), 'r',label='Teoría')
plt.xlabel(r'$\kappa$', fontsize=15)
plt.ylabel('Dimensión fractal', fontsize=12)
plt.grid()
plt.legend()
print('slope = ', lin.slope, '+/-',2*lin.stderr)
print('intercept = ', lin.intercept, '+/-',2*lin.intercept_stderr)
print(f'R^2 = {lin.rvalue**2}')
plt.savefig('k_df')
plt.show()
