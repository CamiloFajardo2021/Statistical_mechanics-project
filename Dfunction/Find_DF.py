import numpy as np
import matplotlib.pyplot as plt
import numexpr as ne
import pandas as pd
import copy as cp
from typing import Type, Union 

Real = Union[np.float16, np.float32, np.float64, float]
Complex = Union[np.complex64, np.complex128] 

def g(zi:Complex, z:Complex) -> Type[Complex]:
    sqrtz = ne.evaluate("sqrt((z - zi.real) ** 2 + zi.imag ** 2)")
    if sqrtz.imag<0:
        return ne.evaluate("zi.real - sqrtz")
    else:
        return ne.evaluate("zi.real + sqrtz")


def var(u:Type[np.array], size:int) -> Type[np.array]:
    var_u = np.zeros(size)
    for ii in range(size):
        u2 = np.dot(u[:ii+1],u[:ii+1])
        mean_u2 = np.mean(u2)
        mean_u = np.mean(u[:ii+1])
        mean_u_2 = np.dot(mean_u, mean_u)
        var_u[ii] = mean_u2-mean_u_2
    return var_u


# Se lee el archivo
Df = pd.read_csv('data2.csv')
Df.columns = ['N', 'zReal', 'zImag']
# Puntos iniciales de z
z = Df['zReal'].to_numpy() +1j*Df['zImag'].to_numpy()
# Arreglo que contiene las transformaciones de z
z_trans = cp.deepcopy(z)

# Arreglo que contiene el tiempo
t = np.zeros(len(z))
# Arreglo que contiene la varianza de la fncion directora 
var_u = np.zeros(len(t))

# Se realiza la trasnformación 
for ii in range(len(z)):
   zi = cp.deepcopy(z_trans[ii])
   t[ii] = (zi.imag**2)/4
   for jj in range(ii, len(z)):
       z_trans[jj]=g(zi, z_trans[jj])

# Se determinan los tiempos
t = np.cumsum(t)
#Se calcula la varianza para cada tiempo
varU = var(z_trans.real, len(t))

#Se imprimen archivos
data = np.column_stack([t, varU])
datafile_path = "data_k_2_25.txt"
np.savetxt(datafile_path , data)

