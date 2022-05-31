import numpy as np
import matplotlib.pyplot as plt

def plot_data(x, y):
    plt.plot(x,y)
    plt.grid()
    plt.show()

#xdata = np.zeros(10)
#xdata = np.linspace(0, 1, 10)
ydata = np.linspace(0, 100, 100)
xdata = np.random.normal(0, 1, 100)
#ydata = np.random.normal(0, 0.5, 100)


n_datos = len(xdata)

deltax = np.zeros(n_datos - 1)
deltay = np.zeros(n_datos - 1)
angles = np.zeros(n_datos - 1)

for ii in range (n_datos - 1):
    deltax[ii] = xdata[ii+1] - xdata[ii]
    deltay[ii] = ydata[ii+1] - ydata[ii]
    if deltax[ii] > 0:
        angles[ii] = np.pi/2 - np.arctan(deltay[ii]/deltax[ii])
    else:
        angles[ii] = np.arctan(deltay[ii]/-deltax[ii]) - np.pi/2
    
print (angles)

print (np.average(angles))

plot_data(xdata, ydata)
