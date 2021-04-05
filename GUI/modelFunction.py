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

    mrna = y[0]
    protein = y[1]

    dmrna = c1*u - c2*mrna
    dprotein = c3*mrna - c4*protein

    return np.array([dmrna, dprotein])