#function created to solve the system of differential equations 
import numpy as np 
from numba import jit 

@jit 
def ODEmodel(y,t,hog,vPars): 
    if t > len(hog): 
        t = len(hog)-1 
    u = hog[int(t)]
    c1 = vPars[0]
    c2 = vPars[1]
    c3 = vPars[2]
    c4 = vPars[3]

    dx1 = y[0]
    dx2 = y[1]
    dx3 = y[2]
    dx4 = y[3]
    dx5 = y[4]

    ddx1 = c1*u - c2*dx1
    ddx2 = c3*dx1 - c4*dx2
    ddx3 = 2*c1*dx1*u + c1*u + c2*dx1 - 2*c2*dx3
    ddx4 = c1*dx2*u - c2*dx4 + c3*dx3 - c4*dx4
    ddx5 = c3*dx1 + 2*c3*dx4 + c4*dx2 - 2*c4*dx5

    return np.array([ddx1, ddx2, ddx3, ddx4, ddx5])