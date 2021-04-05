# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 09:45:58 2020

@author: Bry

#module containing several functions created to simulate and modelate the gen
expression of cell populations

FUNCTIONS:
    1) HOGexpr: Computes the input system in pulses shape.
    
    2) simbODE: Determines simbolic ODEs that describe system output.
        
    3) solveODE: Solves deterministically simbolic expressions of ODEs.
    
    4) gillAl: Solves stochastically simbolic expressions of ODEs using 
        Gillespie Algorithm.
        
    5) cGill: Solves stochastically simbolic expressions of ODEs using 
        Gillespie Algorithm. Main loop is compiled in C.
        
    6) simbMoments: Determines simbolic ODEs that describe the moments of 
        the system up to second order.
        
    7) MLLmeasure: Calculates minus log likelihood measure using observations
        and simulation. 
    
    8) solve2Moment: Solves simbolic expressions of the system moments.
    
    9) kldmeasure: Calculates kld measure from observation and simulation.
        
    10) modelDefiner: creates a .py file that contains the differential 
        equation system.
    
    11) solveODEpy: uses the .py function created by modelDefiner ans computes
        its numerical deterministic solution of the differential equation 
        system.
    
    12) model2MDefiner: Create a .py file that contains the moments of the 
        system
    
    13) solveMODEpy: Solve deterministically raw moments of the system.
    
    14) solve2M. Compute first and second moment of the observed species.
"""

#libraries required in this module
import importlib
import numpy as np
import sympy as sym
from scipy.integrate import odeint
import modelFunction
#import C_Gill as cg
#import matplotlib.pyplot as plt

#1)
########################### pulse_expr ######################################
#inputs: pulse start vector(list), pulse durantion vector(list),
    #duration of simulation
def pulse_expr(ton, tdur, tend):
    ###### VALVE #######
    #time vector
    cell_time = np.linspace(0,int(tend),int(tend) + 1)
    
    #create profile of valve shock
    ### SHOCK VALVE ###
    #check pulse starts are greater than 0
    #checks pulse starts are smaller than simulation duration
    vton = np.where((ton >= 0) & (ton < tend))
    
    ton = ton[vton]
    tdur = tdur[vton]
    
    #Check that pulse duration does not overlap with the next pulse
    for i in range(len(ton) - 1):
        if ton[i]+tdur[i] > ton[i+1]:
            print("Pulses overlap")
    
    uValve = np.zeros((len(cell_time)))

    for i in range(len(ton)):    
        uValve[int(ton[i]):int(ton[i]) + int(tdur[i])] = 1.0

    return uValve, cell_time
#end pulse_expr
##########################################################################
#2)
########################### simbODE ######################################
#inputs: molecular species, reactives matrix, products matrix, parameters 
    #names, input name, index of reaction affected by input

def simbODE(species, reactants, products, pars, inputN=None, indexU=1):
    #molecular species
    species = sym.var(species)
    
    #stoichiometric matrix
    V = products - reactants
    
    #reaction network dimentions
    Sn, Rm = V.shape
    
    #kinetic parameters
    pars = sym.var(pars)
    
    #pre-propensity function
    if inputN == None:
        aPro = pars.copy()
        nameVar = species.copy()
    
    else:
        #system input
        inputN = [inputN]
        uh = sym.var(inputN)
        
        aPro = pars.copy()
        aPro[indexU-1] = aPro[indexU-1]*uh[0]  #include input to propensity
                                                                    #vector
        #variable names
        nameVar = [uh[0]]
        for i in range(0, Sn):
            nameVar.append(species[i])
    #end if
    
    #propensity function vector
    for r in range(0,Rm):
        for s in range(0, Sn):
            #determine propensity vector expressions
            for a in range(0, reactants[s,r]):
                aPro[r] *= species[s]
            #end for a
        #end for s
    #end for r

    #ODEs system
    odeX = []
    for s in range(0,Sn):
        temp = 0

        #determine ODEs system
        for r in range(0, Rm):
            temp += V[s,r]*aPro[r]
        #end for r
    
        #set of ODEs
        odeX.append(temp)
    #end for s
    
    #store variable names in a dictionary
    names = {
            "species":species,
            "pars":np.array(pars),
            "nameVar":np.array(nameVar)
            }

    #return simbolic ODEs and variable names
    return np.array(odeX), names
#end simbODE
    
########################################################################################
#3)    
############################    solveODE   #############################################
#inputs: parameter values, regressors(simbolic ODEs, variable names, input 
    #stimulus, time vector, initial concentrations)

def solveODE(Vpars, inputs):
    #parameters names
    namePars = inputs["pars"]
    
    #simbolic ODEs
    exprs = inputs["ODEs"]
    
    #replace parameter values in simbolic expressions
    odePars = []
    for expr in exprs:
        for i in range(0, len(namePars)):
            expr = expr.subs(namePars[i], Vpars[i])
        #end for i
        odePars.append(expr)
    #end for expr
    
    #ODE model
    def modelODE(z,t,hog):
        #append values to evaluate system
        tempZ = np.array([hog])
        tempZ = np.append(tempZ, z)
        
        #check amount of values to evaluate
        if len(nameVar) != len(tempZ):
            tempZ = tempZ[1:]
        
        #evaluate expressions
        evalue = np.array(exp(tempZ))

        return evalue
    #end modelODE
    
    #variable names
    nameVar = inputs["nameVar"]   
    exp = sym.lambdify([nameVar], odePars, "numpy")
    
    #initial concentrations
    sp0 = inputs["species0"]
    
    #time vector
    tog = inputs["Vtime"]
    
    #input system
    hog = inputs["inpU"]
    
    #array to store solution
    valuesSp = np.zeros((len(sp0), len(tog)))
    
    #include initial concentrations
    valuesSp[:,0] = sp0
    
    #solve ODEs system
    for t in range(1, int(tog[-1]) + 1):
        #time interval
        tspan = np.array([tog[t-1], tog[t]])
        
        #integrate ODEs in time interval
        z = odeint(modelODE,sp0,tspan,args=(hog[t],))
        
        #store computed values
        valuesSp[:,t] = z[1,:]
        
        #change initial contidion for next iteration
        sp0 = z[1]
    #end for t
    
    return valuesSp
    #end solveODE

############################################################################
#4)
#################   gillAL ################################################
#inputs: parameters, regressors(reactives matrix, products matrix, input 
    #stimulus, time vector, initial concentrations, index of reaction affected
    #by input)

def gillAl(parsValues, inputs):    
    #stoichiometric reactions
    reactants = inputs["matrizR"]
    products = inputs["matrizP"]

    #reaction affected by input
    indexU = inputs["idxR"]
    if indexU == None:
        indexU = 0
        inputHog = np.ones((len(inputs["inpU"])))
    else:
        indexU = indexU - 1
        #system input
        inputHog = inputs["inpU"]
    #################################################
    #include fake reaction 
    cfake = np.zeros((reactants.shape[0],1))
    
    reactants = np.hstack((reactants, cfake))
    products = np.hstack((products, cfake))
    
    rfake = np.zeros((1,reactants.shape[1]))
    
    reactants = np.vstack((reactants, rfake))
    products = np.vstack((products, rfake))
    products[-1,-1] = 2.0
    
    reactants = reactants.astype(int)
    products = products.astype(int)
    ###################################################  
    #initial concentrations
    A0 = inputs["species0"] 
    
    #time vector
    tog = inputs["Vtime"]
    
    #create parameters copy and include fake reaction and species
    knu = np.copy(parsValues)
    knu = np.append(parsValues, 1)
    A0 = np.append(A0, 0)
    
    #####################################################################
    #paramters copy that changes during simulation
    real_knu = np.copy(knu)
    
    #take number of species and parameters
    Nchs = A0.shape[0]
    q = knu.shape[0]
    
    #initial conditions of simulation
    m = 0
    A = np.zeros((len(A0), 1))
    A[:,0] = A0
    t = np.zeros(1)
    
    while t[m] <= tog[-1]:
        #temporal variable to store propensity obtained from reactions and
        #species
        alpha = np.zeros(q)
        
        #STEP 1
        #for each parameter calculate propensity taken each previous reaction
        
        #compute av=hv*cv. Where cv are reaction parameters. hv is the number 
        #different molecular reactant in a reaction Ru
        for i in range(0,q):
            if i == indexU and indexU != None:
                approx_t = np.floor(t[m]).astype(int)# + 1
                real_knu[i] = knu[i]*inputHog[approx_t]
            #end if i
        
            alpha[i] = real_knu[i]
    
            for j in range(0,Nchs):
                for s in range(0,reactants[j,i]):
                    ss = 1
                    alpha[i] = alpha[i]*(A[j,m]-ss+1)
                #end for s
            #end for j
            
        #compute a0 = sum(av)          
        alpha0 = np.sum(alpha)
        #end for i
        
        #STEP 2
        #generate two random numbers. First one is use to calculate tau value
        #tau = (1/a0)*ln(1/r1). tau represent instant of time which a reaction
        #takes place
        r = np.random.rand(2)

        alpha0_1 = 1/alpha0
        tau = alpha0_1*(np.log(1/r[0]))
        
        #r2 is used to compute reaction type that takes place for the next 
        #iteration in the instant time tau. ii = u. u is the index of the 
        # --> sum(avi/ao) > r2        
        ii = 0
        sumalpha = alpha0_1*alpha[0]
    
        while r[1] >= sumalpha:
            ii += 1
            sumalpha += alpha0_1*alpha[ii]
        #end while r
        
        #STEP 3
        #adjust time and concentration values 
        #increment one m value
        Atemp = A[:,m] - reactants[:,ii] + products[:,ii] 
        A = np.concatenate((A, np.array([Atemp]).T), axis=1)
        t = np.append(t, (t[m] + tau))
        
        #uptade m
        m += 1
        #while loop is repeadted until a final condition is achieved
    #end while t
    
    #insert final value
    if t[-1] > tog[-1]:
        t[-1] = tog[-1]
    
    #remove values under 0 and fake species
    A[A < 0] = 0
    Apost = np.copy(A[:-1,:])
      
    #new array to store species values
    Npl = 2*(len(t) - 2) + 2
    Apl = np.zeros((Nchs - 1, Npl))
    
    
    #insert repeated values
    for i in range(0,len(t) - 1):#-1
        Apl[:,2*i] = Apost[:,i]
        Apl[:,2*i+1] = Apost[:,i]
    #end for i
    
    #create a new array for time and assign it start and end values
    tpl = np.zeros(Npl) 
    tpl[0] = t[0]
    tpl[Npl-1] = t[len(t) - 1]
    
    #insert repeated values through the whole vector
    for j in range(1,len(t) - 1):
        tpl[2*j-1] = t[j]
        tpl[2*j] = t[j] 
    #end for j
    
    #remove repeated values
    tpl2, ind = np.unique(tpl, return_index=True)
    Apl2 = Apl[:,ind]
    
    #new array to store species profiles
    Apl2x = np.zeros((len(tog),Nchs - 1))

    for inter in range(0,Nchs - 1):
        Apl2x[:,inter] = np.interp(tog, tpl2, Apl2[inter,:].T)
    #end for inter
    
    #remove values under 0
    Apl2x[Apl2x<0] = 0
    cellExpr = Apl2x.T
    
    return cellExpr
#end for cell
##############################################################################
#5)
######################### cGill ##############################################
#inputs: parameters, regressors(reactives matrix, products matrix, input 
    #stimulus, time vector, initial concentrations, index of reaction affected
    #by input)

def cGill(parsValues, inputs):
    #stoichiometric reactions
    reactants = inputs["matrizR"]
    products = inputs["matrizP"]

    #reaction affected by the input stimulus
    indexU = inputs["idxR"]
    if indexU == None:
        indexU = 0
        inputHog = np.ones((len(inputs["inpU"])))
    else:
        indexU = indexU - 1
        #input system
        inputHog = inputs["inpU"]
    #################################################
    #include fake reaction to reactives and products matrices
    cfake = np.zeros((reactants.shape[0],1))
    
    reactants = np.hstack((reactants, cfake))
    products = np.hstack((products, cfake))
    
    rfake = np.zeros((1,reactants.shape[1]))
    
    reactants = np.vstack((reactants, rfake))
    products = np.vstack((products, rfake))
    products[-1,-1] = 2.0
    
    reactants = reactants.astype(int)
    products = products.astype(int)
    ###################################################  
    #initial concentrations
    A0 = inputs["species0"] 
    
    #time vector
    tog = inputs["Vtime"]
    
    #parameters copy containing fake reaction parameter and append fake species
    #initial concentration
    knu = np.copy(parsValues)
    knu = np.append(parsValues, 1)
    A0 = np.append(A0, 0)
    
    #####################################################################
    #parameters copy that changes during simulation
    real_knu = np.copy(knu)
    
    #number of species and kinetica paramers
    Nchs = A0.shape[0]
    q = knu.shape[0]
    
    #preallocation array
    Nalloc = 1000000
    
    A = np.zeros((Nchs, Nalloc))
    t = np.zeros(Nalloc)
    
    #initial conditions to perform simulation
    m = 0
    A[:,0] = A0
    t[0] = 0
    
    #duration of experiment
    tend = tog[-1]
    
    #compiled while loop to solve equations
    Aexpr = cg.gillLoop(A, reactants, products, knu, real_knu, inputHog, tog, 
                        tend, t, m, q, Nchs, indexU)
    
    return Aexpr
##############################################################################
#6)    
###################### simbMoments ###########################################
#inputs: molecular species names, reactives matrix, products matrix, 
    #kinetic parameters names, input stimulus name, reaction affected by input

def simbMoments(species, reactants, products, pars, inputN=None, indexU=1):
    #molecular species
    species = sym.var(species)
    
    #stoichiometric matrix
    V = products - reactants
    
    #reaction network dimentions
    Sn, Rm = V.shape
    
    #kinetica parameters
    pars = sym.var(pars)
    
    #pre-función de propensidad
    #pre-propensity function
    if inputN == None:
        aPro = pars.copy()
    
        nameVar = species.copy()
    
    else:
        #system input
        inputN = [inputN]
        uh = sym.var(inputN)
        
        aPro = pars.copy()
        aPro[indexU-1] = aPro[indexU-1]*uh[0]  #includes input stimulus
        
        #variable names
        nameVar = [uh[0]]
        for i in range(0, Sn):
            nameVar.append(species[i])
    #end if
    
    #propensity vector function
    for r in range(0,Rm):
        for s in range(0, Sn):
            #computes propensity vector
            for a in range(0, reactants[s,r]):
                aPro[r] *= species[s]
            #end for a
        #end for s
    #end for r

    #ODEs system
    odeX = []
    for s in range(0,Sn):
        temp = 0
        
        #determines first order ODEs
        for r in range(0, Rm):
            temp += V[s,r]*aPro[r]
        #end for r
        
        #set of ODEs
        odeX.append(temp)
    #end for s
    
    #empty lists to store ODEs system
    ode2m = []
    name2m = []
    nameODE = species
    odeTotal = odeX
    for s1 in range(0, Sn):
        for s2 in range(0, Sn):
            #determine second order expressions
            temp = 0
            for r in range(0,Rm):
                temp += (V[s1,r]*aPro[r]*species[s2] + V[s2,r]*aPro[r]*species[s1] \
                         + V[s1,r]*V[s2,r]*aPro[r])
            #end for r
            
            #set of ODEs of second order
            if temp not in ode2m:
                ode2m.append(temp)
                odeTotal.append(temp)
            #end if temp
            
            #variable names of second order
            if (species[s1]*species[s2]) not in name2m:
                name2m.append(species[s1]*species[s2])
                nameODE.append(species[s1]*species[s2])
            #end if species s1*s2
        #end for s2
    #end for s1
    
    #replace of variable names
    dxName = []
    dxODE = []
    for exp in odeTotal:
        #each name is changed for a nickname 'dx#'
        for j in range(0,len(nameODE)):
            name = sym.var('dx' + str(len(nameODE)-j))
            exp = exp.subs(nameODE[len(nameODE)-j-1],name)
            
            #save varialbes nicknames
            if name not in dxName:
                dxName.append(name)
            #end if name
        #end for j
        dxODE.append(exp)
    #end for exp
    dxName.reverse()
    
    #variable names
    if inputN == None:
        nameVar = dxName.copy()
    else:
        nameVar = [uh[0]]
        for i in range(0, len(dxName)):
            nameVar.append(dxName[i])
        #end for i
    nameVar = np.array(nameVar)
    dxODE = np.array(dxODE)

    #store variable names in a dictionary
    names = {
            "species":species,
            "spODE":odeTotal,
            "pars":np.array(pars),
            "nameVar":np.array(nameVar)
            }
    #return simbolic ODEs and variable names
    return dxODE, names
#end simbMoments

#############################################################################
#7)    
################  MLLmeasure ################################################
#inputs: observations, simulation, variaility of simulation
def MLLmeasure(Obs, fSim, hSim):
    hSim = np.absolute(hSim)    #avoid negative-zero log errors
    mll = np.sum(np.sum(0.5*((Obs - fSim)/hSim)**2, axis=1)) + np.sum(np.sum(np.log(hSim), axis=1))
    return mll

###############################################################################
#8)
#################   solve2Moment ##############################################
#Inputs: kinetic parameters, noise parameters, regressors (ODEs, variable 
    #names,input stimulus, time vectos, initial conditions
def solve2Moment(parsValues, errV, regressor2):
    #computes raw moments
    moments = solveODE(parsValues, regressor2)

    #noise parameters
    a = errV[0]
    b = errV[1]
    
    #takes first moment corresponding to the output mean
    regressorPre = regressor2["regressor"]
    MBmean = moments[len(regressorPre["species"]) - 1] #1st moment
    
    MBvar = moments[regressor2["idx2M"]] - MBmean**2 #variance 2nd raw moment
    #sd2M = np.sqrt(MBvar) #standard deviation 2nd raw moment
    MBsd_noisy = np.sqrt((1 + b**2)*MBvar + (a + b*MBmean)**2) #standard
    #deviation of second moment of the observed species
    
    return MBmean, MBsd_noisy
#############################################################################
#9)    
###################### KLDmeasure ###########################################
#inputs: mean and standard deviation observations and simulation
    
def KLDmeasure(uObs, sdObs, uM, sdM):
    #avoid negative-zero log errors
    sdM = np.absolute(sdM)
    sdObs = np.absolute(sdObs)
    sdObs[sdObs == 0] = 0.01
    
    kld = np.log(sdM/sdObs) + ((sdObs**2 + (uObs-uM)**2)/(2*sdM**2)) - 1/2
    return np.mean(kld)

###############################################################################
#10)    
##################### modelDefiner ############################################
#creates a .py file from the input data. This file contains the
#differential equations system. modelDefiner takes as input the sympy objects
    #corresponding to the species, differential equations and kinectic paramers
    #names
def modelDefiner(species, odeX, pars, modeD = 0):
    #creates differential equations in string format
    dtVar = ['d' + str(i) for i in species]
    dtODE = [dtVar[i] + ' = ' + str(odeX[i]) for i in range(len(odeX))]
    var0 = [str(species[i]) + ' = y[' + str(i) + ']' for i in range(len(species))]
    dtPars = [str(pars[i]) + ' = vPars[' + str(i) + ']' for i in range(len(pars))]
    
    strVar = ', '.join(dtVar)
    ########################################################
    #creates model
    fname = 'modelFunction.py'
    
    #funcion head
    if modeD == 0:
        fStart = """#function created to solve the system of differential equations \nimport numpy as np \nfrom numba import jit \n\n@jit \ndef ODEmodel(y,t,hog,vPars): \n    if t > len(hog): \n        t = len(hog)-1 \n    u = hog[int(t)]"""
    else:
        fStart = """#function created to solve the system of differential equations \nimport numpy as np \ndef ODEmodel(y,t,hog,vPars): \n    if t > len(hog): \n        t = len(hog)-1 \n    u = hog[int(t)]"""
        
    #function end
    fEnd = '\n\n    return np.array([' + strVar + '])'
    
    #creates .py file that contains the ODE function of the system
    with open(fname, 'w+') as f:
        #creates head of the function definition
        f.write(fStart)
        #adds parameters
        for i in range(len(dtPars)):
            f.write('\n    ' + dtPars[i])
        f.write('\n')
        #adds initial conditions
        for i in range(len(var0)):
            f.write('\n    ' + var0[i])
        f.write('\n')
        #adds differential equations
        for i in range(len(dtODE)):
            f.write('\n    ' + dtODE[i])
        #creates function return
        f.write(fEnd)
    
    #imports model as a global object
    import modelFunction
    global modelFunction
    
    #update module
    try: 
        importlib.reload(modelFunction)
    except:
        pass
##############################################################################
#11)
######################### solveODEpy ##########################################
#input: parameters values, regressors(initial concentrations, time vector,
    #input stimulus)
def solveODEpy(Vpars, inputs):
    #initial concentrations
    sp0 = inputs["species0"]
    
    #time vector
    tog = inputs["Vtime"]
    
    #system input
    hog = inputs["inpU"]
    
    #array to store solution
    valuesSp = np.zeros((len(sp0), len(tog)))
    
    #apppends initial concentrations
    valuesSp[:,0] = sp0
    
    #solves differential equations system
    sol = odeint(modelFunction.ODEmodel, sp0, tog, args=(hog, Vpars))
    sol = interp_nan(sol.transpose())
    
    if np.mean(sol[-1,:]) == 0:
        for t in range(1, int(tog[-1]) + 1):
            #time interval
            tspan = np.array([tog[t-1], tog[t]])
            
            #computes solution in each time interval
            z = odeint(modelFunction.ODEmodel,sp0,tspan,args=(hog, Vpars))
            
            #stores computed values
            valuesSp[:,t] = z[1,:]
            
            #initial concentration for next time interval
            sp0 = z[1]
        sol = valuesSp

    return sol

###############################################################################
#12)
################################# model2MDefiner ##############################
#inputs: species, differential equations, parameters names
def model2MDefiner(species, odeX, pars, modeD = 0):
    #creates differential equations in string format
    dtVar = ['d' + str(i) for i in species]
    dtODE = [dtVar[i] + ' = ' + str(odeX[i]) for i in range(len(odeX))]
    var0 = [str(species[i]) + ' = y[' + str(i) + ']' for i in range(len(species))]
    dtPars = [str(pars[i]) + ' = vPars[' + str(i) + ']' for i in range(len(pars))]
    
    strVar = ', '.join(dtVar)
    ########################################################
    #creates model
    fname = 'model2moment.py'
    
    #funcion head
    if modeD == 0:
        fStart = """#function created to solve the system of differential equations \nimport numpy as np \nfrom numba import jit \n\n@jit \ndef ODEmodel(y,t,hog,vPars): \n    if t > len(hog): \n        t = len(hog)-1 \n    u = hog[int(t)]"""
    else:
        fStart = """#function created to solve the system of differential equations \nimport numpy as np \ndef ODEmodel(y,t,hog,vPars): \n    if t > len(hog): \n        t = len(hog)-1 \n    u = hog[int(t)]"""
    
    #function end
    fEnd = '\n\n    return np.array([' + strVar + '])'
    
    #creates .py file that contains the ODE function of the system
    with open(fname, 'w+') as f:
        #creates head of the function definition
        f.write(fStart)
        #adds parameters
        for i in range(len(dtPars)):
            f.write('\n    ' + dtPars[i])
        f.write('\n')
        #adds initial conditions
        for i in range(len(var0)):
            f.write('\n    ' + var0[i])
        f.write('\n')
        #adds differential equations
        for i in range(len(dtODE)):
            f.write('\n    ' + dtODE[i])
        #creates function return
        f.write(fEnd)
    
    #import created model as a global object
    import model2moment
    global model2moment
    
    #update module
    try: 
        importlib.reload(model2moment)
    except:
        pass
    
#############################################################################
#13)
################### solve2MODEpy ############################################
#inputs: kinetica parameters, regressors (initial concentrations, time vector,
    # input stimulus)

def solve2MODEpy(Vpars, inputs):
    #initial concentrations
    sp0 = inputs["species0"]
    
    #time vector
    tog = inputs["Vtime"]
    
    #system input
    hog = inputs["inpU"]
    
    #array to store solution
    valuesSp = np.zeros((len(sp0), len(tog)))
    
    #apppends initial concentrations
    valuesSp[:,0] = sp0
    
    #solves differential equations system
    sol = odeint(model2moment.ODEmodel, sp0, tog, args=(hog, Vpars))
    sol = interp_nan(sol.transpose())
    
    if np.mean(sol[-1,:]) == 0:
        for t in range(1, int(tog[-1]) + 1):
            #time interval
            tspan = np.array([tog[t-1], tog[t]])
            
            #computes solution in each time interval
            z = odeint(model2moment.ODEmodel,sp0,tspan,args=(hog, Vpars))
            
            #stores computed values
            valuesSp[:,t] = z[1,:]
            
            #initial concentration for next time interval
            sp0 = z[1]
        sol = valuesSp
    
    return sol

##############################################################################
#14)
######################### solve2M ############################################
#inputs: kinetic parameters, noise parameters, regressors (input stimulus, 
    #initial concentration, time vector)
def solve2M(pars, noiseP, regressor):
    #noise parameters
    a = noiseP[0]
    b = noiseP[1]
    
    #solve raw moments
    moments = solve2MODEpy(pars, regressor)
    
    idx1M = len(regressor["regressor"]["species"]) - 1 
    
    #1st moment of the observed species
    uM = moments[idx1M,:]
    #2nd moment of the observed species
    sgM = np.sqrt(moments[-1,:] - moments[idx1M,:]**2)
    sgMn = np.sqrt((1 + b**2)*(sgM**2) + (a + b*uM)**2)
    return uM, sgMn
#############################################################################
    #15)
    #################### interp_nan #########################################
#interpolate nan values in a numpy array
def interp_nan(arrNan):
    def nan_helper(y):
        """Helper to handle indices and logical indices of NaNs.
    
        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """
        return np.isnan(y), lambda z: z.nonzero()[0]

    for i in range(len(arrNan)):
        nans, x = nan_helper(arrNan[i,:])
        arrNan[i,nans] = np.interp(x(nans), x(~nans), arrNan[i,~nans])
        
    return arrNan
