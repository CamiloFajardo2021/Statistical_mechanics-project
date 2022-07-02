import matplotlib.pyplot as plt
import numpy as np
import numexpr as ne
import copy as cp
import pandas as pd
from typing import Type, Union 
from matplotlib.animation import FFMpegWriter, PillowWriter

np.random.seed(14322)
Real = Union[np.float16, np.float32, np.float64, float]
Complex = Union[np.complex64, np.complex128]

def main():
    t = np.linspace(0, 10, 1000)
    for i in range(0,9):
        kappa = float(i)
        u = driving_function(len(t), t[1]-t[0], kappa)
        z = SLE(t, u, kappa)
        #zeros = np.zeros(len(t))

        plt.plot(z.real, abs(z.imag), '-')
        #plt.grid()
        plt.savefig(f"SLE_evolution-final_kappa={kappa: 0.1f}.png")

        #Saving the data in csv
        z_data = {'z real': z.real,
        'z imag': abs(z.imag)}
        z_data = pd.DataFrame(z_data)
        z_data.to_csv(f'final_curve_values_kappa={kappa: 0.1f}.csv')


########################### SECONDARY FUNCTIONS ###########################
def driving_function(size:int, dt:Real, k:Real) -> Type[np.array]:
    u = np.random.standard_normal(size)
    u = np.cumsum(np.sqrt(dt*k)*u)
    return u

def vslit(w:Complex, dt:Real, u:Real) -> Type[Complex]:
    return ne.evaluate("1j * sqrt(4 * dt - (w-u) ** 2) + u")

def SLE(t:Type[np.array], u:Type[np.array], kappa:Real) -> Type[np.array]:
    
    #Animation parameters
    metadata = dict(title = "SLE_evolution", artist = "Group 5 - Statistical Mechanics")
    writer = FFMpegWriter(fps = 30, metadata = metadata)
    fig = plt.figure(figsize = (7,7))
    plt.xlim(-1, 4.5)
    plt.ylim(0, 6.5)
    l, = plt.plot([],[], "-k", linewidth = 1)

    #SLE evolution
    nsteps = len(t)
    z = np.zeros(len(t), dtype=np.complex_)
    z = cp.deepcopy(u+0.00000001j)
    
    with writer.saving(fig, f"SLE_evolution-final_kappa={kappa: 0.1f}.mp4", dpi = 100):
        for step in range(nsteps-1, 0, -1):
            dt = t[step]-t[step-1]
            z[step:] = vslit(z[step:], dt, float(u[step]))
            l.set_data(z.real, abs(z.imag))

            writer.grab_frame()

    return z

if __name__ == "__main__":
    main()


