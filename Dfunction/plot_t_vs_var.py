import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def fun (x, a):
    return a*x


t1, varU1 = np.loadtxt('./data/data_k_1.txt', unpack=True)
t2, varU2 = np.loadtxt('./data/data_k_1_25.txt', unpack=True)
t3, varU3 = np.loadtxt('./data/data_k_1_5.txt', unpack=True)
t4, varU4 = np.loadtxt('./data/data_k_1_75.txt', unpack=True)
t5, varU5 = np.loadtxt('./data/data_k_2.txt', unpack=True)
t6, varU6 = np.loadtxt('./data/data_k_2_25.txt', unpack=True)


popt1, pcov1 = curve_fit(fun, t1, varU1/10)
popt2, pcov2 = curve_fit(fun, t2, varU2/10)
popt3, pcov3 = curve_fit(fun, t3, varU3/10)
popt4, pcov4 = curve_fit(fun, t4, varU4/10)
popt5, pcov5 = curve_fit(fun, t5, varU5/10)
popt6, pcov6 = curve_fit(fun, t6, varU6/10)

plt.grid()
plt.ylabel(r'$\langle U_{t}^{2}\rangle-\langle U_{t}\rangle^{2}$', fontsize=12)
plt.xlabel('t', fontsize=14)
'''
plt.plot(t1, varU1/10, 'ko', markersize=3)
plt.plot(t1, t1, 'r--', label='Teórica')
plt.title(r'$k=1.0$', fontsize=15)
plt.legend()
plt.savefig('t_vs_var_k_1')
'''
'''
plt.plot(t2, varU2/10, 'bo', markersize=3)
plt.plot(t2, 1.25*t2, 'r--', label='Teórica')
plt.title(r'$k=1.25$', fontsize=15)
plt.legend()
plt.savefig('t_vs_var_k_1_25')
'''
'''
plt.plot(t3, varU3/10, 'go', markersize=3)
plt.plot(t3, 1.5*t3, 'r--', label='Teórica')
plt.title(r'$k=1.5$', fontsize=15)
plt.legend()
plt.savefig('t_vs_var_k_1_5')
'''
'''
plt.plot(t4, varU4/10, 'ro', markersize=3)
plt.plot(t4, 1.75*t4, 'k--', label='Teórica')
plt.title(r'$k=1.75$', fontsize=15)
plt.legend()
plt.savefig('t_vs_var_k_1_75')
'''
'''
plt.plot(t5, varU5/10, 'mo', markersize=3)
plt.plot(t5, 2.0*t5, 'r--', label='Teórica')
plt.title(r'$k=2.0$', fontsize=15)
plt.legend()
plt.savefig('t_vs_var_k_2')
'''
'''
plt.plot(t6, varU6/10, 'yo', markersize=3)
plt.plot(t6, 2.25*t6, 'r--', label='Teórica')
plt.title(r'$k=2.25$', fontsize=15)
plt.legend()
plt.savefig('t_vs_var_k_2_25')
#plt.show()
'''
print('Slope', '\t','Ds')
print(popt1, 3*np.sqrt(np.diag(pcov1)))
print(popt2, 3*np.sqrt(np.diag(pcov2)))
print(popt3, 3*np.sqrt(np.diag(pcov3)))
print(popt4, 3*np.sqrt(np.diag(pcov4)))
print(popt5, 3*np.sqrt(np.diag(pcov5)))
print(popt6, 3*np.sqrt(np.diag(pcov6)))
