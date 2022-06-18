import numpy as np
import matplotlib.pyplot as plt
import cmath
import math

# Numerical simulation of chordal SLE_k over the interval [0,T]
k = 3.0
T = 0.1

# Tolerence used by the variable step size control
tol = 0.0125

# Maximum and minimum step sizes allowed in the simulation
max_stepsize = T/2048.0
min_stepsize = T/(2.0**33.0)

sqrt_k = math.sqrt(k)

# Lists containing the discretized SLE trace and Brownian path.
sle_path = [0.0+0.0j]
brownian_path = []

# Computes the SLE trace associated with a constant driver.
# This is achieved by solving z' = - 2h/z at t = 0.5.
def horizontal_trace(z0, h):
    zt = cmath.sqrt(z0**2.0 - 2.0*h)
    
    if zt.imag < 0:
        zt = -zt
        
    return zt


# Computes the SLE trace associated with a vertical driver.
# This is achieved by solving z' = sqrt_k*increment at t = 1.
def vertical_trace(z0, increment):

    return z0 + sqrt_k*increment


# Computes the SLE trace using a step of the Ninomiya-Victoir scheme.
# This is equivalent to solving the backward Loewner equation driven
# by a piecewise linear path that has horizontal and vertical pieces.
def ninomiya_victoir(z0, h, brownian_increment):
    
    return horizontal_trace(vertical_trace(horizontal_trace(z0, h),
                                           brownian_increment), h)


# Propagates the numerical SLE with a variable step-size so that
#
# |SLE_{i+1} - SLE_{i}| < tol    for all neighbouring discretization points.
#
# This type of step size control was proposed for chordal SLE in the paper:
# Tom Kennedy, Numerical Computations for the Schramm-Loewner Evolution,
# Journal of Statistical Physics, 137:839, 2009.
def sle_step(time_increment, brownian_increment):
    l = len(brownian_path)
    
    # Approximate the next SLE point using the Ninomiya-Victoir scheme.
    # SLE traces are obtained by solving Loewner's differential equation
    # backwards in time starting from zero (hence the use of reversed()).
    zt = ninomiya_victoir(0.0+0.0j, time_increment, brownian_increment)
    
    for m in reversed(range(l)):
        zt = ninomiya_victoir(zt, brownian_path[m][0], brownian_path[m][1])
        
    # Check the propsed SLE point is sufficiently close to the previous point
    if ((abs(zt - sle_path[l-1]) < tol) or (time_increment <= min_stepsize))  \
    and (time_increment <= max_stepsize):
        
        # Update the numerical SLE trace and Brownian path        
        brownian_path.append([time_increment, brownian_increment])    
        sle_path.append(zt)

    else:
        # Generate the midpoint of the Brownian path over
        # the time interval using the Brownian bridge.        
        bridge_midpoint = 0.5*np.random.normal(0.0, math.sqrt(time_increment))

        half_time_increment = 0.5*time_increment
        half_increment = 0.5*brownian_increment
        
        # Approximate the SLE trace over the two subintervals
        zt = sle_step(half_time_increment, half_increment + bridge_midpoint)
        zt = sle_step(half_time_increment, half_increment - bridge_midpoint)
        
    return zt


# Compute the numerical SLE trace over the interval [0,T] using sle_step
zt = sle_step(T, np.random.normal(0.0, math.sqrt(T)));

# Plot the numerical SLE trace using matplotlib.pyplot
no_of_points = len(sle_path)

xvalues = []
yvalues = []

for i in range(no_of_points):
    xvalues.append(sle_path[i].real)
    yvalues.append(sle_path[i].imag)

    """
plt.figure(0, dpi=300)
plt.plot(xvalues, yvalues, linewidth=0.5)
plt.title('SLE with k = ' + str(round(k, 2)))
plt.axis('scaled')
plt.show()
"""

# Find the index of the Ly parameter
def Ly_index(ydata):
    maxpos = ydata.index(max(ydata))
    return maxpos


# Find the Ly parameter
def Ly_value(ydata):
    return max(ydata)

    
# Find the alpha_i between 3 vectors
def angle_alpha_i(x1,x2,x3,y1,y2,y3):
    v1 = np.array([x2-x1,y2-y1])
    v2 = np.array([x3-x2,y3-y2])
    mag1 = np.sqrt(np.dot(v1, v1))
    mag2 = np.sqrt(np.dot(v2, v2))
    cose = abs(np.dot(v1,v2)/(mag1*mag2))
    #theta = math.degrees(math.acos(cose))
    theta = math.acos(cose)
    
    def slope():
        m1 = (y2-y1)/(x2-x1)
        m2 = (y3-y2)/(x3-x2)
        return m1, m2
    
    if (y2 > y1) and (x2 > x3):
        theta = -theta 
    
    elif (y1 > y2) and (x3 > x2):
        theta = -theta
    
    elif (y1 == y2):
        if (y2 > y3) and (x1 > x2):
            theta = -theta
        if (y3 > y2) and (x2 > x1):
            theta = -theta
    else:
        m1, m2 = slope()
        if m2 > m1:
            theta = -theta


    return theta



# Find the winding angle for all datapoints
def winding_angle(xdata, ydata):
    n = len(xdata)
    angles = np.zeros(n-1)
    angles[0] = 0
    #angles[0] = math.atan((ydata[1]-ydata[0])/(xdata[1]-xdata[0]))
    for ii in range (1, n-1):
        alpha_i = angle_alpha_i(xdata[ii-1], xdata[ii], xdata[ii+1],
                               ydata[ii-1], ydata[ii], ydata[ii+1])
        angles[ii] = angles[ii-1] + alpha_i
        #angles[ii] = alpha_i
    
    return angles
    


# Find the average of the squares of a vector
# <theta**2> 
def average_squares(vector):
    n = len(vector)
    suma = 0
    for ii in range (n):
        suma += (vector[ii]**2)
    return suma/(n)


# Trick to not simulate over and over
# Sebas dice que no sirve xd
def portion(cte, xdata, ydata):
    n = len(xdata)
    xdata = xdata[0: int(cte*n): 1]
    ydata = ydata[0: int(cte*n): 1]
    return xdata, ydata

    
print('Number of steps = ' + str(no_of_points-1))
print("\n \nLy and <theta**2> are")


angles = winding_angle(xvalues, yvalues)
average_theta2 = average_squares(angles)
Ly = Ly_value(yvalues)
print (Ly, average_theta2)

