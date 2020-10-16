# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 09:45:58 2020

@author: Bry

Modulo con funciones para la simulación de la expresión genica en poblaciones
celulares

FUNCIONES:
    1) HOGexpr: Calcula la expresión del hog como entrada del sistema.
    
    2) simbODE: Calculo simbolico de las ODEs que describen la respuesta del
        sistema.
        
    3) solveODE: Soluciona las expresiones simbolicas de las ODEs.
    
    4) simGill: Calculo de las expresion estocastica de cada una de las 
        especies de forma especifica para cada celula perteneciente a una
        población
        
    5) cGill: Calculo de las expresion estocastica de cada una de las 
        especies de forma especifica para cada celula perteneciente a una
        población. Función con ciclo while compilado en Cython.
        
    6) simbMoments: Calculo simbolico del sistemas de ecuaciones diferenciales
        que definen la respuesta de las especies moleculares y de sus 
        respectivos momentos de segundo orden.
        
    7) MLLmeasure: Calcula la medida de verosimilitud logaritmica negativa 
        entre las observaciones y la simulación. Usado para modelo de Célula 
        Promedio. Por lo tanto, la simulación hace referencia a la simulación de una
        misma célula.
    
    8) solve2Moment: Soluciona las expresiones simbolicas de los momentos de 
        segundo ordern del sistema.
    
    9) kldmeasure: Calcula la medida kld a partir de las observaciones y las 
        simulaciones. Toma datos estadisticos de media y desvicación estandar.
"""

#librerias
import numpy as np
#import symengine as sym
#https://github.com/symengine/symengine.py
import sympy as sym
from scipy.integrate import odeint
import C_Gill as cg
#import matplotlib.pyplot as plt

########################### HOGexpr ######################################
#entradas: tiempos de inicio pulso (list), tiempo de duración de pulso(list),
    #tiempo de duracion de la simulación
def HOGexpr(ton, tdur, tend):
    ###### VALVE #######
    #Vector de tiempo del experimento
    cell_time = np.linspace(0,int(tend),int(tend) + 1)
    
    #Creación del perfil del Shock de la valvula
    ### SHOCK VALVE ###
    #Verificación de que el inicio de los pulsos sea mayor que 0
    #Verificación de que el inicio de los pulsos sea menor al instante en que 
    #termina el experimento
    vton = np.where((ton >= 0) & (ton < tend))
    
    ton = ton[vton]
    tdur = tdur[vton]
    
    #Se comprueba que la duración de los pulsos no se sobreponga con el pulso 
    #siguiente
    for i in range(len(ton) - 1):
        if ton[i]+tdur[i] > ton[i+1]:
            print("Los pulsos de la valvula se sobreponen")
    
    #Valores iniciales de tiempo y del shock        
    t0 = np.array([0], float)
    u0 = np.array([0], float)
    
    #Valor alto del shock en los instantes de inicio de los pulsos
    uon = ton*0 + 1
    
    #Valor alto del shock en los instantes de terminación de los pulsos
    toff = ton + tdur
    uoff = 0*toff + 1
    
    #En el mismmo instante de inicio del shock, el estado del shock al principio
    #estaba en bajo
    tpre = ton
    upre = 0*tpre
    
    #En el instante de terminación del shock, el estado del shock pasa a bajo
    tpos = toff
    upos = 0*tpos
    
    #Estado final del shock 
    if toff[-1] > tend:
        uend = np.array([uoff[-1]])
    else:
        uend = np.array([0], float)
    
    #Concatenación de los instantes de tiempo y sus respectivos estados durante el
    #experimento
    tValve = np.concatenate((t0, tpre, ton, toff, tpos, tend))
    uValve = np.concatenate((u0, upre, uon, uoff, upos, uend))
    
    #Se sacan los indices correspodientes al orden ascendentes de los instantes de
    #tiempo
    ind = np.argsort(tValve)
    
    #Se toman los indices para ordenar los vectores de tiempo y de sus 
    #correspondientes estados del shock de la valvula
    tValve = tValve[ind]
    uValve = uValve[ind]
    
    #Verificación y eliminación de que algun instante de tiempo no supere el valor 
    #de tiempo de fin del experimento
    nonv = np.where(tValve > tend[0])
    
    tValve = np.delete(tValve, nonv)
    uValve = np.delete(uValve, nonv)
    
    #Arreglo que contiene el comportamiento de la valvula
    t_u_Valve = np.vstack((tValve,uValve))
    
#    print(t_u_Valve)
#    plt.figure()
#    plt.plot(tValve,uValve)
#    plt.show()
    ######  CHAMBER ##############
    #Delay en el inicio del perfil del shock en la camara
    ton = ton + 3
    
    #Evolución del perfil del shock en la camara
    #Valores en alto en los instantes de inicio del shock en la camara
    uon = ton*0 + 1
    
    #valores en alto hasta los instantes de fin del shock. Se le resta un instante
    #en el cual la evolución del perfil en la camara desciende de forma exponencial
    toff = ton + tdur - 1
    uoff = 0*toff + 1
    
    #El perfil del shock en la camara en un instante previo al inicio del shock
    #toma valor de bajo
    tpre = ton - 1
    upre = 0*tpre
    
    #La decaida del perfil de la camara evoluciona en cuatro instantes de tiempo
    tpos1 = toff + 1
    upos1 = tpos1*0 + 0.75
    
    tpos2 = toff + 2
    upos2 = tpos2*0 + 0.5
    
    tpos3 = toff + 3
    upos3 = tpos3*0 + 0.25
    
    tpos4 = toff + 4
    upos4 = tpos4*0
    
    #Correccion del error que cambia los valores de inicio del siguiente por los de
    #parte posterior del shock actual
    
    #Vectores temporales con los instantes de tiempo y valores de estado de la camara
    #posteriores
    temp0 = np.concatenate((tpos1, tpos2, tpos3, tpos4))
    utemp0 = np.concatenate((upos1, upos2, upos3, upos4))
    
    #Toma los instantes de tiempo posteriores correspondientes a cada shock y los compara con 
    #inicio del shock siguiente. De este modo, se conserva en alto si se llega a dar una sobreposición
    #de los perfiles en la camara
    for i in range(0, len(uon) - 1):
        tempt = temp0[i::(len(uon))]
        
        for j in range(0,4):
            if tempt[j] > ton[i+1]:
                if j == 0:
                    upos1[i] = 1
                elif j == 1:
                    upos2[i] = 1
                elif j == 2:
                    upos3[i] = 1
                elif j == 3:
                    upos4[i] = 1
        #end for j
    #end for i 
    
    #Correccion del cambio realizado por el estado previo del perfil de un shock
    #Evita  que la parte posterior del perfil actual se ponga en bajo por la parte
    #previa del perfil siguiente
    for i in range(0, len(uon) - 1):
        tempt = temp0[i::(len(uon))]
        tempu = utemp0[i::(len(uon))]
        
        for j in range(0,4):
            if tpre[i+1] == tempt[j]:
                upre[i+1] = tempu[j]
        #end for j
    #end for i
    
    #Valor final del estado del perfil
    if toff[-1] > tend:
        uend = np.array([uoff[-1]])
    else:
        uend = np.array([0], float)
    
    #Concateanción de los instantes de tiempo y los estados correspondientes
    #de la evolución del perfil de los shock en al camara
    tChamber = np.concatenate((t0, tpre, ton, toff, tpos1, tpos2, tpos3, tpos4, tend))
    uChamber = np.concatenate((u0, upre, uon, uoff, upos1, upos2, upos3, upos4, uend))
    
    #Se ordenan los instantes de tiempo de forma ascendente y se toman los indices
    # de los instantes sin repetición
    uniq, ind = np.unique(tChamber, return_index=True)
    
    #A partir de los indices se toman los valores del estado de la camara y sus 
    #respectivos instantes de tiempo
    tChamber = tChamber[ind]
    uChamber = uChamber[ind]
    
    #se verifica y se eliminan los instantes de tiempo de la camara que superen el instante de
    #fin del experimento
    nonc = np.where(tChamber > tend[0])
    
    tChamber = np.delete(tChamber, nonc)
    uChamber = np.delete(uChamber, nonc)
    
    #Arreglo que contiene el comportamiento del perfil del shock en la camara
    t_u_Chamber = np.vstack((tChamber,uChamber))
    
#    plt.figure()
#    plt.subplot(2,1,1)
#    plt.plot(tValve, uValve, 'r-',label='Valve(t)')
#    plt.legend(loc='best')
#    
#    plt.subplot(2,1,2)
#    plt.plot(tChamber, uChamber, 'g-',label='Chamber(t)')
#    plt.legend(loc='best')
#    plt.show()
    
    ####### EXPRESION DEL HOG  ##########
    t_u_Valve = np.interp(cell_time, t_u_Valve[0], t_u_Valve[1])
    u = np.interp(cell_time, tChamber, uChamber)
    
    #modelo del HOG
    def modelHOG(z,t,u):
        hog = z[0]
        dydt = 0.3968*u - 0.9225*hog
        return dydt
    
    #valor inicial de HOG
    hog0 = 0.0
    
    #creación de vector que contendra los valores de HOG
    hog = np.zeros(len(cell_time))
    hog[0] = hog0
    
    #solución de la ODE del sistema
    for i in range(1,int(tend) + 1):
        tspan = [cell_time[i-1], cell_time[i]]
        
        z = odeint(modelHOG,hog0,tspan,args=(u[i],))
        
        hog[i] = z[1]
        hog0 = z[1]
    #end for i
    
    ################################
#    plt.figure()
#    plt.subplot(3,1,1)
#    plt.plot(tValve, uValve, 'r-',label='Valve(t)')
#    plt.legend(loc='best')
#    
#    plt.subplot(3,1,2)
#    plt.plot(tChamber, uChamber, 'g-',label='Chamber(t)')
#    plt.legend(loc='best')
#    
#    plt.subplot(3,1,3)
#    plt.plot(cell_time, hog, 'b-',label='HOG(t)')
#    plt.legend(loc='best')
#    plt.show()
    
    #Diccionario con perfiles de la valvula y la camara
    profiles = {
            "t_u_Valve":t_u_Valve,
            "t_u_Chamber":t_u_Chamber
            }
    
    #regresa la expresión del hog, vector de tiempo del experimento, y 
    #perfiles previos a la expresion
    return hog, cell_time, profiles
#end HOGexpr
    
#TEST
#CREACIOIN DEL HOG COMO ENTRADA AL SISTEMA

##Se define el instante de tiempo en el cual termina el experimento
#tend = np.array([100], float)
#
##Arreglo con los instantes de tiempo de inicio de los shocks y sus respectivas
##duraciones
#ton = np.array([1], float)
#tdur = np.array([3], float)
#
#inputHog, tog, perfiles = HOGexpr(ton, tdur, tend)
##########################################################################

########################### simbODE ######################################
#Entradas: Nombres especies moleculares(list, string), matrix de reactivos, 
        #matrix de productos, #Nombres parametros cineticos (list, string), 
        #Opcional: Nombre entrada del sistema (string)
        #Reacción que afecta la entrada (index)

def simbODE(species, reactants, products, pars, inputN=None, indexU=1):
    #Especies moleculares
    species = sym.var(species)
    
    #Matriz estequiometrica
    V = products - reactants
    
    #dimensiones  de la red de reacciones
    Sn, Rm = V.shape
    
    #parametros cineticos
    pars = sym.var(pars)
    
    #pre-función de propensidad
    if inputN == None:
        aPro = pars.copy()
        nameVar = species.copy()
    
    else:
        #Entrada del sistema
        inputN = [inputN]
        uh = sym.var(inputN)
        
        aPro = pars.copy()
        aPro[indexU-1] = aPro[indexU-1]*uh[0]  #añade la entrada a la reacción
        #correspondiente
        
        #nombres de las variables
        nameVar = [uh[0]]
        for i in range(0, Sn):
            nameVar.append(species[i])
    #end if
    
    #vector de función de propensidad
    for r in range(0,Rm):
        for s in range(0, Sn):
            #calculo de las expresiones del vector de propensidad
            for a in range(0, reactants[s,r]):
                aPro[r] *= species[s]
            #end for a
        #end for s
    #end for r

    #Sistema de ODEs
    odeX = []
    for s in range(0,Sn):
        temp = 0
        #calculo de la ODE
        for r in range(0, Rm):
            temp += V[s,r]*aPro[r]
        #end for r
        #conjunto de ODEs del sistema
        odeX.append(temp)
    #end for s
    
    #diccionario con los nombres de las variables
    names = {
            "species":species,
            "pars":np.array(pars),
            "nameVar":np.array(nameVar)
            }
    #Regresa arreglo con las ODEs simbolicas, y nombres de las variables
    return np.array(odeX), names
#end simbODE
    
##TEST
#especies = ['X1', 'X2']
#
#reactivos = np.array([[0, 1, 1, 0],[0, 0, 0, 1]])
#productos = np.array([[1, 0, 1, 0],[0, 0, 1, 0]])
#
#parametros = ['c1', 'c2', 'c3', 'c4']
#
#entrada = 'U'
#
##ecuaciones, variables = simbODE(especies, reactivos, productos, parametros)
#ecuaciones, variables = simbODE(especies, reactivos, productos, parametros, inputN=entrada)
#print(variables)
#print(ecuaciones)
########################################################################################
    
############################    solveODE   #############################################
#Entradas: valores de los parametros, Entradas(ODEs, nombreVariales, inpU, time,
    #condiciones iniciales)

def solveODE(Vpars, inputs):
    
    #nombre de las variables parametros
    namePars = inputs["pars"]
    
    #ODEs simbolicas
    exprs = inputs["ODEs"]
    
    #reemplazo de los valores de los parametros en las expresiones simbolicas.
    odePars = []
    for expr in exprs:
        for i in range(0, len(namePars)):
            expr = expr.subs(namePars[i], Vpars[i])
        #end for i
        odePars.append(expr)
    #end for expr
    
    #Modelo ODE
    def modelODE(z,t,hog):
        #adjunta datos a evaluar
        tempZ = np.array([hog])
        tempZ = np.append(tempZ, z)
        
        #verifica el número de datos a evaluar            
        if len(nameVar) != len(tempZ):
            tempZ = tempZ[1:]
        
        #Evalua las expresiones            
        evalue = np.array(exp(tempZ))

        return evalue
    #end modelODE
    
    #nombre de las variables
    nameVar = inputs["nameVar"]   
    exp = sym.lambdify([nameVar], odePars, "numpy")
    
    #condiciones iniciales de las especies
    sp0 = inputs["species0"]
    
    #vector de tiempo
    tog = inputs["Vtime"]
    
    #entrada del sistema
    hog = inputs["inpU"]
    
    #arreglo para almacenar la solución
    valuesSp = np.zeros((len(sp0), len(tog)))
    
    #registra condiciones iniciales
    valuesSp[:,0] = sp0
    
    #solución del sistema de ODEs
    for t in range(1, int(tog[-1]) + 1):
        #intervalo de tiempo
        tspan = np.array([tog[t-1], tog[t]])
        
        #solución de las ODEs en cada intervalo de tiempo
        z = odeint(modelODE,sp0,tspan,args=(hog[t],))
        
        #almacena los valores calculados
        valuesSp[:,t] = z[1,:]
        
        #condicion inicial del siguiente intervalo
        sp0 = z[1]
    #end for t
    
    return valuesSp
    #end solveODE

#regresa los valores de concentracion de cada una de las especies.
    
##parametros cineticos
#parsValues = [4.0, 0.010, 1.0, 0.006]
#
##condiciones iniciales de las especies
#sp0 = np.zeros(len(especies))
#
#
#Allins = {
#        "ODEs":ecuaciones,
#        "inpU":inputHog,
#        "Vtime":tog,
#        "species0":sp0
#        }
#
#Allins.update(variables)
#
#exprEspecies = solveODE(parsValues, Allins)
############################################################################

#################   simGill ################################################
#entradas: parametros, reactivos, productos, entradas del
    #sistema (time, input, condiciones iniciales, indice de reaccion que 
    #afecta la entrada)

def gillAl(parsValues, inputs):
    
    #Matrices estequiometricas de reacciones
    reactants = inputs["matrizR"]
    products = inputs["matrizP"]

    #reaccion que afecta la entrada
    indexU = inputs["idxR"]
    if indexU == None:
        indexU = 0
        inputHog = np.ones((len(inputs["inpU"])))
    else:
        indexU = indexU - 1
        #entrada del sistema
        inputHog = inputs["inpU"]
    #################################################
    #adicion de la especie y la reacción falsa
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
    #Condiciones iniciales
    A0 = inputs["species0"] 
    
    #vector de tiempo
    tog = inputs["Vtime"]
    
    #Se crea otra variable con los parametros. A esta variable se le añade un parametro
    #adicional para la especie falsa. De igual modo, se añade la concentración inicial
    #de dicha especie
    knu = np.copy(parsValues)
    knu = np.append(parsValues, 1)
    A0 = np.append(A0, 0)
    
    #####################################################################
    #Esta variale contiene igualmente los valores de los parametros. Sin embargo, 
    #en este caso el valor del parametro que es afectado por la entrada
    #cambiara en el tiempo dependiendo del valor de entrada
    real_knu = np.copy(knu)
    
    #Se toman los valores del numero de especies y el numero de parametros  
    Nchs = A0.shape[0]
    q = knu.shape[0]
    
    #Condiciones iniciales para la simulación estocastica    
    m = 0
    A = np.zeros((len(A0), 1))
    A[:,0] = A0
    t = np.zeros(1)
    
    while t[m] <= tog[-1]:
        #Se crea una variable temporal para almacenar los calculos de las propensiones
        #obtenidas  de las reacciones y las especies
        alpha = np.zeros(q)
        
        #STEP 1
        
        #Para cada parametro se calculo la propensión con respecto a cada reacción
        #previo, para parametros que cambio con respecto a la entrada se hace
        #la respectiva modificación
        
        #Se calculo av=hv*cv. Donde cv es el parametro de la reacción. hv es el numero
        # de distintos reactantes moleculares en un reacción Ru
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
            
        #Calculo de a0 = sum(av)          
        alpha0 = np.sum(alpha)
        #end for i
        
        #STEP 2
        #Se generan dos numeros aleatorios. Con el primer r1 se calcula el valor de tau
        #tau = (1/a0)*ln(1/r1). tau representa el instante de tiempo en el cual se va 
        #dar la siguiente reaccion
        r = np.random.rand(2)

        alpha0_1 = 1/alpha0
        tau = alpha0_1*(np.log(1/r[0]))
        
        #con r2 se calculo el tipo de reaccion que se va a dar en el instante de tiempo
        #con tau. ii = u. u va a indicar el tipo de reaccion que se dara. Su valor se 
        #determina a partir de la condición --> sum(avi/ao) > r2        
        ii = 0
        sumalpha = alpha0_1*alpha[0]
    
        while r[1] >= sumalpha:
            ii += 1
            sumalpha += alpha0_1*alpha[ii]
        #end while r
        
        #STEP 3
        #se ajustan los valores de tiempo y de las concentraciones de las especies.
        #Se incrementa el valor de m
        Atemp = A[:,m] - reactants[:,ii] + products[:,ii] 
        A = np.concatenate((A, np.array([Atemp]).T), axis=1)
        t = np.append(t, (t[m] + tau))
        
        #Actualiza m
        m += 1
        #El while se repetira hasta que el vector de tiempo supere el valor maximo de tiempo
        #establecido para la simulación
    #end while t
    
    #Se añade el valor final del experimento al vector de tiempo    
    if t[-1] > tog[-1]:
        t[-1] = tog[-1]
    
    #Se eliminan los valores inferiores a 0 y se eliminan los valores de 
    #la especie falsa    
    A[A < 0] = 0
    Apost = np.copy(A[:-1,:])
    
    #Se crea un arreglo con estas dimensiones para almacenar los valores de las 
    #concentraciones de las especies    
    Npl = 2*(len(t) - 2) + 2
    Apl = np.zeros((Nchs - 1, Npl))
    
    #Se añaden valores repetidos en instantes de tiempo consecutivos
    for i in range(0,len(t) - 1):#-1
        Apl[:,2*i] = Apost[:,i]
        Apl[:,2*i+1] = Apost[:,i]
    #end for i
    
    #Se crea un nuevo arreglo para el tiempo, y se asignan los valores de inicio y
    #fin del tiempo simulado 
    tpl = np.zeros(Npl) 
    tpl[0] = t[0]
    tpl[Npl-1] = t[len(t) - 1]
    
    #Se añaden instantes de tiempo repetivos consecutivamente
    for j in range(1,len(t) - 1):
        tpl[2*j-1] = t[j]
        tpl[2*j] = t[j] 
    #end for j
    
    #Se eliminan los valores repetidos y se toman los indices originales de los 
    #valores que quedaron
    tpl2, ind = np.unique(tpl, return_index=True)
    Apl2 = Apl[:,ind]
    
    #Se crea un nuevo arreglo para las concentraciones de las especies.
    #en este caso almacenara la interpolación de los arreglos A y t, en los 
    #instantes de tiempo reales del experimento
    Apl2x = np.zeros((len(tog),Nchs - 1))

    for inter in range(0,Nchs - 1):
        Apl2x[:,inter] = np.interp(tog, tpl2, Apl2[inter,:].T)
    #end for inter
    
    #Se eliminan valors por debajo de 0
    Apl2x[Apl2x<0] = 0
    cellExpr = Apl2x.T
    
    return cellExpr
#end for cell
###############################################################    

######################### cGill ###############################
#entradas: parametros, reactivos, productos, entradas del
#sistema (time, input, condiciones iniciales, indice de reaccion que 
#afecta la entrada)
#toma los mismos datos de entrada que la función anterior. La diferencia
    # es que esta función tiene el ciclo while compilado en cython

def cGill(parsValues, inputs):
        
    #Matrices estequiometricas de reacciones
    reactants = inputs["matrizR"]
    products = inputs["matrizP"]

    #reaccion que afecta la entrada
    indexU = inputs["idxR"]
    if indexU == None:
        indexU = 0
        inputHog = np.ones((len(inputs["inpU"])))
    else:
        indexU = indexU - 1
        #entrada del sistema
        inputHog = inputs["inpU"]
    #################################################
    #adicion de la especie y la reacción falsa
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
    #Condiciones iniciales
    A0 = inputs["species0"] 
    
    #vector de tiempo
    tog = inputs["Vtime"]
    
    #Se crea otra variable con los parametros. A esta variable se le añade un parametro
    #adicional para la especie falsa. De igual modo, se añade la concentración inicial
    #de dicha especie
    knu = np.copy(parsValues)
    knu = np.append(parsValues, 1)
    A0 = np.append(A0, 0)
    
    #####################################################################
    #Esta variale contiene igualmente los valores de los parametros. Sin embargo, 
    #en este caso el valor del parametro que es afectado por la entrada
    #cambiara en el tiempo dependiendo del valor de entrada
    real_knu = np.copy(knu)
    
    #Se toman los valores del numero de especies y el numero de parametros  
    Nchs = A0.shape[0]
    q = knu.shape[0]
    
    #prealocación de arreglos
    Nalloc = 1000000
    
    A = np.zeros((Nchs, Nalloc))
    t = np.zeros(Nalloc)
    
    #Condiciones iniciales para la simulación estocastica    
    m = 0
    A[:,0] = A0
    t[0] = 0
    
    #tiempo de fin del experimento
    tend = tog[-1]
    
    #ciclo while compilado para calcular la expresión
    Aexpr = cg.gillLoop(A, reactants, products, knu, real_knu, inputHog, tog, 
                        tend, t, m, q, Nchs, indexU)
    
    return Aexpr
##############################################################
    
###################### simbMoments ###########################
#Entradas: Nombres especies moleculares(list, string), matrix de reactivos, 
        #matrix de productos, #Nombres parametros cineticos (list, string), 
        #Opcional: Nombre entrada del sistema (string)
        #Reacción que afecta la entrada (index)

def simbMoments(species, reactants, products, pars, inputN=None, indexU=1):
    #Especies moleculares
    species = sym.var(species)
    
    #Matriz estequiometrica
    V = products - reactants
    
    #dimensiones  de la red de reacciones
    Sn, Rm = V.shape
    
    #parametros cineticos
    pars = sym.var(pars)
    
    #pre-función de propensidad
    if inputN == None:
        aPro = pars.copy()
    
        nameVar = species.copy()
    
    else:
        #Entrada del sistema
        inputN = [inputN]
        uh = sym.var(inputN)
        
        aPro = pars.copy()
        aPro[indexU-1] = aPro[indexU-1]*uh[0]  #añade la entrada a la reacción
        #correspondiente
        
        #nombres de las variables
        nameVar = [uh[0]]
        for i in range(0, Sn):
            nameVar.append(species[i])
    #end if
    
    #vector de función de propensidad
    for r in range(0,Rm):
        for s in range(0, Sn):
            #calculo de las expresiones del vector de propensidad
            for a in range(0, reactants[s,r]):
                aPro[r] *= species[s]
            #end for a
        #end for s
    #end for r

    #Sistema de ODEs
    odeX = []
    for s in range(0,Sn):
        temp = 0
        #calculo de la ODE
        for r in range(0, Rm):
            temp += V[s,r]*aPro[r]
        #end for r
        #conjunto de ODEs del sistema
        odeX.append(temp)
    #end for s
    
    #sistema de ODEs de 2do orden
    ode2m = []
    name2m = []
    nameODE = species
    odeTotal = odeX
    for s1 in range(0, Sn):
        for s2 in range(0, Sn):
            #calculo de la expresión de 2 orden
            temp = 0
            for r in range(0,Rm):
                temp += (V[s1,r]*aPro[r]*species[s2] + V[s2,r]*aPro[r]*species[s1] \
                         + V[s1,r]*V[s2,r]*aPro[r])
            #end for r
            #Conjunto de ODEs de momento de 2do orden
            if temp not in ode2m:
                ode2m.append(temp)
                odeTotal.append(temp)
            #end if temp
            #Nombres de las variables de 2do orden
            if (species[s1]*species[s2]) not in name2m:
                name2m.append(species[s1]*species[s2])
                nameODE.append(species[s1]*species[s2])
            #end if species s1*s2
        #end for s2
    #end for s1
    
    #sustitución de nombres de las variables
    dxName = []
    dxODE = []
    for exp in odeTotal:
        #en cada expresión se busca el nombre de la variable
        #para reemplazarlas por un sobrenombre
        #sobrenombre 'dx#'
        for j in range(0,len(nameODE)):
            name = sym.var('dx' + str(len(nameODE)-j))
            exp = exp.subs(nameODE[len(nameODE)-j-1],name)
            
            #guarda los sobrenombre
            if name not in dxName:
                dxName.append(name)
            #end if name
        #end for j
        dxODE.append(exp)
    #end for exp
    dxName.reverse()
    
    #nombres de las variables
    if inputN == None:
        nameVar = dxName.copy()
    else:
        nameVar = [uh[0]]
        for i in range(0, len(dxName)):
            nameVar.append(dxName[i])
        #end for i
    nameVar = np.array(nameVar)
    dxODE = np.array(dxODE)

    #diccionario con los nombres de las variables
    names = {
            "species":species,
            "spODE":odeTotal,
            "pars":np.array(pars),
            "nameVar":np.array(nameVar)
            }
    #Regresa arreglo con las ODEs simbolicas, y nombres de las variables
    return dxODE, names
#end simbMoments

#############################################################
    
################    minus-log likelihood #####################
    #función que calcula la medida de la  de verosimilitud logaritmica
    #negativa
#Entradas: Observaciones, simulación, variabilidad de simulación
def MLLmeasure(Observaciones, fSim, hSim):    
    MLL1 = np.sum(np.sum(0.5*(((Observaciones - fSim)/hSim)**2), axis=1))
    MLL2 = np.sum(np.sum(np.log(hSim), axis=1))
    MLL = MLL1 + MLL2
    return MLL

################################################################
    
#################   solve2Moment ###############################
#Entradas: parametros cineticos, parametors de ruido, regressor con 
    #(ODEs, nombreVariales, inpU, time, condiciones iniciales)
    
def solve2Moment(parsValues, errV, regressor2):
    #calcula los momentos crudos del sistema
    moments = solveODE(parsValues, regressor2)
    
    #calculo de los momentos del sistema
    #parametros de ruido de la medicion
    a = errV[0]
    b = errV[1]
    
    #toma como 1rt momento la salida con respecto a la especie p. 
    #equivale a la respuesta media del sistema
    regressorPre = regressor2["regressor"]
    MBmean = moments[len(regressorPre["species"]) - 1] #1st moment
    
    MBvar = moments[regressor2["idx2M"]] - MBmean**2 #varianza 2nd raw moment
    #sd2M = np.sqrt(MBvar) #desviación estandar 2nd raw moment
    MBsd_noisy = np.sqrt((1 + b**2)*MBvar + (a + b*MBmean)**2) #desviación estandar
    #del segundo momento de la respuesta
    
    return MBmean, MBsd_noisy
####################################################################
    
###################### kld measure #################################
#calcula la medida de kld a partir de la media y la desviación estandar de 
    #las observaciones y las simulaciones
    
def kldmeasure(CMEmean, CMEsd, meanMC2, sdMC2):
    sdMC2 = np.absolute(sdMC2)
    CMEsd = np.absolute(CMEsd)
    CMEsd[CMEsd == 0] = 0.01
    #calculo de la medida de KLD. su medida relaciona la simulación con las observaciones
    kld = np.log(sdMC2/CMEsd) + ((CMEsd**2) + (CMEmean - meanMC2)**2)/(2*sdMC2**2) - 1/2
    #regresa un vector con la medida de KLD para cada instante de tiempo
    return kld