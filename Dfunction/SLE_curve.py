import matplotlib.pyplot as plt
import numpy as np
import numexpr as ne
import copy as cp
import pandas as pd
from typing import Type, Union
import csv

Real = Union[np.float16, np.float32, np.float64, float];
Complex = Union[np.complex64, np.complex128];

np.random.seed(14322)

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


def var(u:Type[np.array], size:int) -> Type[np.array]:
    var_u = np.zeros(size)
    for ii in range(size):
        u2 = np.dot(u[:ii+1],u[:ii+1])
        mean_u2 = np.mean(u2)
        mean_u = np.mean(u[:ii+1])
        mean_u_2 = np.dot(mean_u, mean_u)
        var_u[ii] = mean_u2 -mean_u_2
    return var_u


tf = 0.1
t = np.linspace(0, tf, 200)
kappa = 2.25
u = driving_function(len(t), t[1]-t[0], kappa)
z = SLE(t, u)
z_data = {'z_real':z.real-z.real[0], 'z_imag':z.imag}
z_data = pd.DataFrame(z_data)
z_data.to_csv('data2.csv')
zeros = np.zeros(len(u[:4]))

#plt.plot(z.real-z.real[0], (z.imag), '-')
plt.plot(t, var(u, len(t))/10, 'o')
plt.plot(t, kappa*t)
#plt.grid()
#plt.plot(u[:4]-u[0], zeros, 'o')
plt.show()

