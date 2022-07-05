import matplotlib.pyplot as plt
import numpy as np
import numexpr as ne
import copy as cp
import pandas as pd
from typing import Type, Union
import csv

Real = Union[np.float16, np.float32, np.float64, float];
Complex = Union[np.complex64, np.complex128];

np.random.seed(14722)

def driving_function(size:int, dt:Real, k:Real) -> Type[np.array]:
    u = np.random.standard_normal(size)
    u = np.cumsum(np.sqrt(dt*k)*u)
    return u

def vslit(w:Complex, dt:Real, u:Real) -> Type[Complex]:
    return ne.evaluate("1j * sqrt(4 * dt - (w-u) ** 2) + u")

def SLE(t:Type[np.array], u:Type[np.array]) -> Type[np.array]:

    nsteps = len(t)
    z = np.zeros(len(t), dtype=np.complex_)
    z = cp.deepcopy(u+0.00000001j)
    for step in range(nsteps-1, 0, -1):
        dt = t[step]-t[step-1]
        z[step:] = vslit(z[step:], dt, (u[step]))

    return z


tf = 1.0
t = np.linspace(0, tf, 100000)
kappa = 1.75
u = driving_function(len(t), t[1]-t[0], kappa)
z = SLE(t, u)
z_data = {'z_real':z.real, 'z_imag':z.imag}
z_data = pd.DataFrame(z_data)
z_data.to_csv(f'final_curve_values_t={tf: 0.1f}.csv')
#zeros = np.zeros(len(t))

plt.plot(z.real, (z.imag), '-')
#plt.grid()
plt.show()

