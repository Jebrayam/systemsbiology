# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 07:45:55 2020

@author: Bry
"""
#import libraries used to programm GUI
import sys, re
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QDesktopWidget, QDialog
from PyQt5.QtGui import QIcon
#package to load the .ui file
from PyQt5 import uic
#libraries for processing data
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import fmin
import sympy as sym
import simsysbio as s2b

#build the dialog window used to define the input model
class inputWindow(QDialog):
    def __init__(self):
        QDialog.__init__(self)
        uic.loadUi("inputWindow.ui", self)
        
        #connect button  to define input model
        self.defineInputButton.clicked.connect(self.defineInput)
        #validate information from inputs
        self.inpModelSp.textChanged.connect(self.valiSInp)
        self.inpModelPar.textChanged.connect(self.valiPInp)
        self.inpModelR.textChanged.connect(self.valiRInp)
        self.inpModelP.textChanged.connect(self.valiPrInp)
    
    #visual validation of inputs
    def valiSInp(self):
        val = re.match('^[a-zA-Z0-9, ]+$', self.inpModelSp.text(), re.I)
        if self.inpModelSp.text() == "":
            self.inpModelSp.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inpModelSp.setStyleSheet("border: 2px solid red;")
        else:
            self.inpModelSp.setStyleSheet("border: 2px solid green;")
    
    def valiPInp(self):
        val = re.match('^[-0-9., ]+$', self.inpModelPar.text(), re.I)
        if self.inpModelPar.text() == "":
            self.inpModelPar.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inpModelPar.setStyleSheet("border: 2px solid red;")
        else:
            self.inpModelPar.setStyleSheet("border: 2px solid green;")
            
    def valiRInp(self):
        val = re.match('^[0-1;\n, ]+$', self.inpModelR.toPlainText(), re.I)
        if self.inpModelR.toPlainText() == "":
            self.inpModelR.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inpModelR.setStyleSheet("border: 2px solid red;")
        else:
            self.inpModelR.setStyleSheet("border: 2px solid green;")

    def valiPrInp(self):
        val = re.match('^[0-1;\n, ]+$', self.inpModelP.toPlainText(), re.I)
        if self.inpModelP.toPlainText() == "":
            self.inpModelP.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inpModelP.setStyleSheet("border: 2px solid red;")
        else:
            self.inpModelP.setStyleSheet("border: 2px solid green;") 
    
    #take data to define input model when "Define" is clicked
    def defineInput(self):
        #take the name of the molecular species
        self.SpStrInp = self.inpModelSp.text()
        self.SpStrInp = self.SpStrInp.split(',')
        self.SpStrInp = [x.strip() for x in self.SpStrInp]
        
        #take parameter values
        self.parStrInp = self.inpModelPar.text()
        self.parStrInp = re.findall('[0-9.0-9]+', self.parStrInp)
        
        #take reagent matrix
        self.reacStrInp = self.inpModelR.toPlainText()
        self.reacStrInp = self.reacStrInp.split(';')
        self.reacStrInp = [x.split(',') for x in self.reacStrInp]
        
        #take product matrix
        self.prodStrInp = self.inpModelP.toPlainText()
        self.prodStrInp = self.prodStrInp.split(';')
        self.prodStrInp = [x.split(',') for x in self.prodStrInp]
        
        try:
            #convert parameters and get the name from each one
            self.valParInp = list(map(float, self.parStrInp))
            self.nparStrInp = [('k'+ str(i+1)) for i in range(len(self.parStrInp))]
            
            #convert data from reagent matrix
            self.valReacInp = list()
            for i in self.reacStrInp:
                self.valReacInp.append(list(map(int, i)))
            self.valReacInp = np.array(self.valReacInp)
            print(self.valReacInp)
            
            #convert data from product matrix
            self.valProdInp = list()
            for i in self.prodStrInp:
                self.valProdInp.append(list(map(int, i)))
            self.valProdInp = np.array(self.valProdInp)
            print(self.valProdInp)
            
            #check input dimentions
            if self.valReacInp.shape != self.valProdInp.shape:
                _mainWindow.outODEs.append("\nThe dimensions of the matrices are not equal.") 
                return
            elif (len(self.valReacInp) != len(self.SpStrInp)) or (len(self.valProdInp) != len(self.SpStrInp)):
                _mainWindow.outODEs.append("\nThe number of rows in the matrices is not equal to the number of species.") 
                return
            elif (self.valReacInp.shape[1] != len(self.valParInp)) or (self.valProdInp.shape[1] != len(self.valParInp)):
                _mainWindow.outODEs.append("\nThe number of columns in the matrices is not equal to the number of parameters.") 
                return
            #input model's name
            self.nameInp = 'U'
            
            #compute the simbolic expression of the input model
            self.ecuacionesInp, self.variablesInp = s2b.simbODE(self.SpStrInp,
                        self.valReacInp, self.valProdInp, self.nparStrInp, 
                        inputN=self.nameInp)
            
            #replace the kinetic parameter values in the model
            self.namePars = self.variablesInp["pars"]
            
            self.ecuacionesInpVal = []
            for expr in self.ecuacionesInp:
                for i in range(0, len(self.namePars)):
                    expr = expr.subs(self.namePars[i], self.valParInp[i])
                #end for i
                self.ecuacionesInpVal.append(expr)
            #end for expr
            
            print(self.ecuacionesInpVal)
            #show simbolic expression of the input model on the window
            _mainWindow.outODEs.append("\nInput Model:") 
            for j in range(len(self.ecuacionesInp)):
                    _mainWindow.outODEs.append('\t' + 'd' + self.SpStrInp[j] + '/dt:  ' + str(self.ecuacionesInpVal[j]))
            
            #save data got from input model in a dictionary
            self.regressorInp = {
                "ODEs": self.ecuacionesInp,
                "species0": np.zeros((len(self.ecuacionesInp))),
                "ParValues": self.valParInp
                }
            self.regressorInp.update(self.variablesInp)
                
        except:
            _mainWindow.outODEs.append("\nError defining input model.")

#main window's constructor
class mainWindow(QMainWindow):
    #constructor method of the class
    def __init__(self):
        #start QMainWindow object
        QMainWindow.__init__(self)
        #load the configuration from .ui file in the object
        uic.loadUi("genexpsim.ui", self)
        #window configuration
        self.UIsetting()
        #set plot windows
        self.plotSetting()
        #connect buttons to methods
        self.buttons()
        #initial state of some variables
        self.idxReaction = None
        #dictionary to save data that plays as regressor
        self.regressor = dict()
        #create window instance of the model input definition
        self.inputModelW = inputWindow()
        #disable inputs that depend on other options
        self.setButton.setEnabled(False)
        self.spinAval.setEnabled(False)
        self.spinBval.setEnabled(False)
        self.covMatrix.setEnabled(False)
        self.changePlotButton.setEnabled(False)
        #initial state of the secondary species subplot
        self.subplotSP = False
        
    def UIsetting(self):
        #set icon
        self.setWindowIcon(QIcon('icon.png'))  
        #set window size
        self.setMinimumSize(self.frameSize().width(), 
                            self.frameSize().height())
        self.setMaximumSize(self.frameSize().width(), 
                            self.frameSize().height())
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
    
    def plotSetting(self):
        self.graphPlot.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot.canvas.axes.grid()
        
        self.graphPlot4.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot4.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot4.canvas.axes.grid()
                
        self.graphPlot2.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot2.canvas.axes.set_ylim([0, 1])
        self.graphPlot2.canvas.axes.grid()
        
        self.graphPlot3.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot3.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot3.canvas.axes.grid()

        self.graphPlot5.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot5.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot5.canvas.axes.grid()

    #bind buttons to corresponding actions
    def buttons(self):
        #connexion to buttons
        self.identiButton.clicked.connect(self.identify)
        self.simButton.clicked.connect(self.simulate)
        self.optButton.clicked.connect(self.inferate)
        self.cancelButton.clicked.connect(self.cancelTool)
        self.changePlotButton.clicked.connect(self.changePlot)
        #enable secondary species plot
        self.radioPlot2.clicked.connect(self.enablePlot2)
        #connexion to set model button
        self.setButton.clicked.connect(self.setInputModel)
        #input system radio button
        self.radioEntrada.clicked.connect(self.enable_entrada)
        #enable inputs of variability sources
        self.extrinsicBox.stateChanged.connect(self.enable_cov)
        self.noiseBox.stateChanged.connect(self.enable_noise)
        #spinbox button to get reaction-input index
        self.Ridx.valueChanged.connect(self.changeIdxR)
        #spinbox button to define experiment duration
        self.spinTdur.valueChanged.connect(self.changeTdur)
        #spinbox buttons to define noise parameters
        self.spinAval.valueChanged.connect(self.changeErrorPar)
        self.spinBval.valueChanged.connect(self.changeErrorPar)
        #spinbox button to define cell numeber
        self.numCells.valueChanged.connect(self.changeCells)
        #combobox to choose input type
        self.inputType.currentIndexChanged.connect(self.changeInput)
        #check input information
        self.inEspecies.textChanged.connect(self.valiS)
        self.inPars.textChanged.connect(self.valiP)
        self.inReact.textChanged.connect(self.valiR)
        self.inProd.textChanged.connect(self.valiPr)
        self.concentSp.textChanged.connect(self.valiConcent)
        self.tOnPulse.textChanged.connect(self.valiOn)
        self.tDurPulse.textChanged.connect(self.valiDur)
        self.covMatrix.textChanged.connect(self.valiMcov)
        self.pars0.textChanged.connect(self.valiPar0)
        #combobox to choose inference model
        self.modelBox.currentIndexChanged.connect(self.changeModel)
    
    #input checks
    def valiMcov(self):
        val = re.match('^[0-9.,;\n ]+$', self.covMatrix.toPlainText(), re.I)
        if self.covMatrix.toPlainText() == "":
            self.covMatrix.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.covMatrix.setStyleSheet("border: 2px solid red;")
        else:
            self.covMatrix.setStyleSheet("border: 2px solid green;")

    def valiDur(self):
        val = re.match('^[0-9, ]+$', self.tDurPulse.text(), re.I)
        if self.tDurPulse.text() == "":
            self.tDurPulse.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.tDurPulse.setStyleSheet("border: 2px solid red;")
        else:
            self.tDurPulse.setStyleSheet("border: 2px solid green;")

    def valiOn(self):
        val = re.match('^[0-9, ]+$', self.tOnPulse.text(), re.I)
        if self.tOnPulse.text() == "":
            self.tOnPulse.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.tOnPulse.setStyleSheet("border: 2px solid red;")
        else:
            self.tOnPulse.setStyleSheet("border: 2px solid green;")

    def valiConcent(self):
        val = re.match('^[0-9., ]+$', self.concentSp.text(), re.I)
        if self.concentSp.text() == "":
            self.concentSp.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.concentSp.setStyleSheet("border: 2px solid red;")
        else:
            self.concentSp.setStyleSheet("border: 2px solid green;")
    
    def valiS(self):
        val = re.match('^[a-zA-Z0-9, ]+$', self.inEspecies.text(), re.I)
        if self.inEspecies.text() == "":
            self.inEspecies.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inEspecies.setStyleSheet("border: 2px solid red;")
        else:
            self.inEspecies.setStyleSheet("border: 2px solid green;")
    
    def valiP(self):
        val = re.match('^[-0-9., ]+$', self.inPars.text(), re.I)
        if self.inPars.text() == "":
            self.inPars.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inPars.setStyleSheet("border: 2px solid red;")
        else:
            self.inPars.setStyleSheet("border: 2px solid green;")
        
    def valiPar0(self):
        val = re.match('^[-0-9., ]+$', self.pars0.text(), re.I)
        if self.pars0.text() == "":
            self.pars0.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.pars0.setStyleSheet("border: 2px solid red;")
        else:
            self.pars0.setStyleSheet("border: 2px solid green;")
            
    def valiR(self):
        val = re.match('^[0-1;\n, ]+$', self.inReact.toPlainText(), re.I)
        if self.inReact.toPlainText() == "":
            self.inReact.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inReact.setStyleSheet("border: 2px solid red;")
        else:
            self.inReact.setStyleSheet("border: 2px solid green;")

    def valiPr(self):
        val = re.match('^[0-1;\n, ]+$', self.inProd.toPlainText(), re.I)
        if self.inProd.toPlainText() == "":
            self.inProd.setStyleSheet("border: 2px solid yellow;")
        elif not val:
            self.inProd.setStyleSheet("border: 2px solid red;")
        else:
            self.inProd.setStyleSheet("border: 2px solid green;") 
    
    #action executed when "cancel" is clicked
    def cancelTool(self):
        self.outODEs.clear()
        self.labeliter.clear()
        self.enableSimOpt()
        self.clearPlot()
    
    #clean all plots
    def clearPlot(self):
        self.graphPlot.canvas.axes.clear()
        self.graphPlot2.canvas.axes.clear()
        self.graphPlot3.canvas.axes.clear()
        self.graphPlot4.canvas.axes.clear()
        self.graphPlot5.canvas.axes.clear()

        self.plotSetting()

        self.graphPlot.canvas.draw()
        self.graphPlot2.canvas.draw()
        self.graphPlot3.canvas.draw()
        self.graphPlot4.canvas.draw()
        self.graphPlot5.canvas.draw()

    #enable simulation and inference section
    def enableSimOpt(self, eSim = False, eOpt = False):
        #enable simulation
        self.simButton.setEnabled(eSim)
        self.tabSim.setEnabled(eSim)
        
        #Enable input configuration
        if self.idxReaction == None:
            self.inputType.setEnabled(False)
            self.tOnPulse.setEnabled(False)
            self.tDurPulse.setEnabled(False)
        else:
            self.inputType.setEnabled(True)
            self.tOnPulse.setEnabled(True)
            self.tDurPulse.setEnabled(True)
        
        #enable inference
        self.optButton.setEnabled(eOpt)
        self.tabOpt.setEnabled(eOpt)
        
    #action executed when "define" is clicked
    def identify(self):
        #inhabilita otros botones o limpia entradas
        #disable other buttons or clean inputs
        self.enableSimOpt()
        self.labeliter.clear()
        self.clearPlot()
        self.outPlot.clear()
        
        #take molecular species names
        self.SpStr = self.inEspecies.text()
        self.SpStr = self.SpStr.split(',')
        self.SpStr = [x.strip() for x in self.SpStr]
        
        #take parameter values
        self.parStr = self.inPars.text()
        self.parStr = re.findall('[0-9.0-9]+', self.parStr)
        
        #take reagent matix values
        self.reacStr = self.inReact.toPlainText()
        self.reacStr = self.reacStr.split(';')
        self.reacStr = [x.split(',') for x in self.reacStr]
        
        #take product matrix values
        self.prodStr = self.inProd.toPlainText()
        self.prodStr = self.prodStr.split(';')
        self.prodStr = [x.split(',') for x in self.prodStr]
                
        try:
            #show data added for definig system
            self.showInfo(self.SpStr, self.parStr)
            
            #convert parameters and get parameter names
            self.valPar = list(map(float, self.parStr))
            self.nparStr = [('c'+ str(i+1)) for i in range(len(self.parStr))]
            print(self.valPar)
            print(self.nparStr)
            print(self.SpStr)
            
            #convert data from reagent matrix
            self.valReac = list()
            for i in self.reacStr:
                self.valReac.append(list(map(int, i)))
            self.valReac = np.array(self.valReac)
            print(self.valReac)
            
            #convert data from product matrix
            self.valProd = list()
            for i in self.prodStr:
                self.valProd.append(list(map(int, i)))
            self.valProd = np.array(self.valProd)
            print(self.valProd)
            print(self.idxReaction)
            
            #check input dimentions
            if self.valReac.shape != self.valProd.shape:
                self.outODEs.append("\nThe dimensions of the matrices are not equal.") 
                return
            elif (len(self.valReac) != len(self.SpStr)) or (len(self.valProd) != len(self.SpStr)):
                self.outODEs.append("\nThe number of rows in the matrices is not equal to the number of species.") 
                return
            elif (self.valReac.shape[1] != len(self.valPar)) or (self.valProd.shape[1] != len(self.valPar)):
                self.outODEs.append("\nThe number of columns in the matrices is not equal to the number of parameters.") 
                return

            #compute simbolic expressions of the system
            #get ODEs in simbolic form together to their respective varible names
            self.nameInp = 'U'
            
            #compute system without input
            if self.idxReaction == None:
                self.ecuaciones, self.variables = s2b.simbODE(self.SpStr,
                    self.valReac, self.valProd,self.nparStr)
                
            #Compute system with input
            else:
                self.ecuaciones, self.variables = s2b.simbODE(self.SpStr,
                    self.valReac, self.valProd, self.nparStr, 
                    inputN=self.nameInp, indexU = self.idxReaction)
            
            #show added data for defining
            self.showInfo(self.SpStr, self.parStr, self.valReac, self.valProd,
                          self.ecuaciones)
            
            #active button: Simulation
            self.enableSimOpt(True)
            
            #store used information in simulation
            self.regressor = {
                    "ODEs": self.ecuaciones,
                    "matrizR": self.valReac,
                    "matrizP": self.valProd,
                    "vPars": self.valPar,
                    "idxR": self.idxReaction
                    }
            self.regressor.update(self.variables)
            
            #set species names in combobox
            self.SpStr.reverse()
            self.outPlot.addItems(self.SpStr)
            self.SpStr.reverse()
            
            #noise measurement messagge
            self.noiseBox.setToolTip("Measurement noise will affect \na the species: " + str(self.SpStr[-1]))
            
        except:
            self.outODEs.append("\nPlease review the input information.")
            
    def showInfo(self, SpStr, parStr, valReac = np.array([]), valProd = np.array([]), ecuaciones = np.array([])):
            
            #show molecular species
            self.outODEs.setText("Molecular Species: " + ', '.join(self.SpStr))
            
            #create parameter names
            if parStr[0] != '':
                nparStr = [('c'+ str(i+1)) for i in range(len(self.parStr))]
                nparStrV = [(nparStr[i] + ' = ' + str(self.parStr[i])) for i in range(len(self.parStr))]
                self.outODEs.append("Kinetic Parameters: " + ', '.join(nparStrV))
            else:
                self.outODEs.append("Kinetic Parameters: ")
            
            #show stoichiometrix matrix
            if valReac.size == 0 or valProd.size == 0:
                return
            else:
                MatrizV = valProd - valReac
                self.outODEs.append("Stoichiometric Matrix: ")
                for i in range(len(MatrizV)):
                    temp = MatrizV[i].tolist()
                    temp = list(map(str, temp))
                    self.outODEs.append('\t\t|  ' + ', '.join(temp) + '  |  ' + SpStr[i])
                print(MatrizV)

            #shoe differential equations
            self.outODEs.append("Differential Equations: ")
            for j in range(len(ecuaciones)):
                self.outODEs.append('\t' + 'd' + SpStr[j] + '/dt =  ' + str(ecuaciones[j]))
                    
    #take reaction-input index if it changes
    def changeIdxR(self):
        self.idxReaction = self.Ridx.value()

    def enable_entrada(self):
        #enable or disable inputs for definig of one system input
        self.parStr = None
        self.parStr = self.inPars.text()
        self.parStr = self.parStr.split(',')
        
        #enable input just if it detects any parameter
        if self.parStr[0] != '':
            if self.radioEntrada.isChecked(): 
                self.Ridx.setEnabled(True)
                self.Ridx.setRange(1, len(self.parStr))
                self.idxReaction = self.Ridx.value()
            else: 
                self.Ridx.setEnabled(False)
                self.idxReaction = None
        else: 
            self.Ridx.setEnabled(False)
            self.idxReaction = None
    
    #action executed when "simulate" is clicked
    def simulate(self):
        #reboot current and previous sections
        self.clearPlot()
        self.changePlotButton.setEnabled(False)
        self.spPlotBox.clear()
        self.spPlotBox2.clear()

        #disable during simulation
        self.enableSimOpt(False)
        
        #take experiment duration value
        self.expTdur = self.spinTdur.value()
        print(self.expTdur)
        
        #number of cells to simulate
        self.nCells = self.numCells.value()
        
        #take initial concentrations values
        self.concenStr = self.concentSp.text()
        self.concenStr = re.findall('[0-9.0-9]+', self.concenStr)
        
        #take input configuration if there is one
        if self.idxReaction != None:
            #take pulse start values
            self.tOnStr = self.tOnPulse.text()
            self.tOnStr = re.findall('[0-9.0-9]+', self.tOnStr)
        
            #take pulse duration values
            self.tDurStr = self.tDurPulse.text()
            self.tDurStr = re.findall('[0-9.0-9]+', self.tDurStr)
        
#        try:
        #convert initial concentration values to float
        self.valConcen = list(map(float, self.concenStr))
        self.regressor["species0"] = np.array(self.valConcen)
    
        #check concentration values to be equal to amount of species
        if len(self.SpStr) != len(self.valConcen):
            self.outODEs.append("\nThe number of initial concentration values is not equal to the number of species.")
            self.enableSimOpt(True)
            return
        print(self.valConcen)
        
        #compute input signal
        if self.idxReaction == None:
            #time and input vectors
            self.expTime = np.linspace(0, self.expTdur, self.expTdur + 1)
            self.inpShock = np.ones(len(self.expTime))
            
        else:
            #convet pulse values to float
            self.valTon = list(map(float, self.tOnStr))
            self.valTdur = list(map(float, self.tDurStr))
            print(self.valTon, self.valTdur)
            
            #check experiment duration
            if not(all((i >= 0 and i < self.expTdur) for i in self.valTon)):
                self.outODEs.append("\nThe start of the pulses is not within the established limits.")
                self.enableSimOpt(True)
                return
            elif not(all((i > 0) for i in self.valTdur)):
                self.outODEs.append("\nThe duration of the pulses must be greater than zero.")
                self.enableSimOpt(True)
                return
            
            #get input type
            self.changeInput()
            print(self.inputT, type(self.inputT))
            
            #compute input system
            self.expTime, self.inpShock = self.computeInput()
        
        #save information for simulating
        self.regressor["Vtime"] = self.expTime
        self.regressor["inpU"] = self.inpShock
        print(self.regressor)
        
        #show information for simulating
        self.showInfosim(self.expTdur, self.concenStr, self.nCells)
        
        #compute system simulation
        self.AllcellExpr = self.computeSim(self.nCells, self.regressor)
        try:
            print(self.AllcellExpr.shape)
        except:
            print(self.AllcellExpr)
            
        #plot system output
        self.idxPlot = self.SpStr.index(self.outPlot.currentText())
        self.plotting(self.expTime, self.inpShock, self.AllcellExpr, self.idxPlot)
        
        #store simulated data
        self.data = {
                "genExpr": self.AllcellExpr,
                "inPulse": self.inpShock,
                "vTime": self.expTime,
                }
        
        #enable button and combobox to change species plot
        self.changePlotButton.setEnabled(True)
    
        #set species names in combobo which change plot
        self.SpStr2 = self.SpStr.copy()
        self.SpStr2.append('All')
        self.spPlotBox.addItems(self.SpStr2)
        self.spPlotBox2.addItems(self.SpStr)
        
        #enable inference section
        self.enableSimOpt(True, True)

#        except:
#            self.enableSimOpt(True)
#            self.outODEs.append("\nPor favor, revisar la información de entrada")
    
    def computeInput(self):
        #compute input
        hog, expTime, perfiles = s2b.HOGexpr(np.array(self.valTon),
                                             np.array(self.valTdur), 
                                             np.array([self.expTdur], float))
        if self.inputT == 'Pulse':
            #take valve values as input
            inpShock = perfiles["t_u_Valve"]
            print(inpShock.shape)
            
        elif self.inputT == 'Model':
            self.inputModelW.regressorInp["inpU"] = perfiles["t_u_Valve"]
            self.inputModelW.regressorInp["Vtime"] = expTime
            print(self.inputModelW.regressorInp)
            print(perfiles["t_u_Valve"].shape, expTime.shape)
            
            #compute system input
            ExprInp = s2b.solveODE(self.inputModelW.regressorInp["ParValues"], self.inputModelW.regressorInp)
            inpShock = ExprInp[-1,:]
            
        return expTime, inpShock
           
    def plotting(self, expTime, inpShock, AllcellExpr, idxPlot): 
        #set if there is a secondary plot
        if self.radioPlot2.isChecked():
            self.graphPlot.resize(891,271)
            self.graphPlot3.resize(891,271)

        else:
            self.graphPlot.resize(891,411)
            self.graphPlot3.resize(891,411)
            
        #show plotting state
        self.labeliter.setText("Plotting...")
        QApplication.processEvents()
        
        #plot input system
        self.graphPlot2.canvas.axes.plot(expTime, inpShock, 'r')
        self.graphPlot2.canvas.axes.fill_between(expTime, inpShock, color='r')
        
        #get species to be plot
        print(self.outPlot.currentText())
        print(idxPlot)
        
        #plot when there is variability
        if self.noiseBox.isChecked() or self.intrinsicBox.isChecked() or self.extrinsicBox.isChecked():
            for i in range(0, self.nCells):
                self.graphPlot.canvas.axes.plot(expTime, AllcellExpr[i,:,idxPlot])
                
                #secondary species plot
                if self.radioPlot2.isChecked() and self.subplotSP == True:
                    self.graphPlot4.canvas.axes.plot(expTime, AllcellExpr[i,:,self.indexPlot2])  
                
                self.labeliter.setText("Plotting... " + str(i+1))
                QApplication.processEvents()
        
            #adjust data to be contain between the 95 % of the output. (2.5-97.5)
            populationlim = np.zeros((AllcellExpr[:,:,idxPlot].shape))
            
            #set percentiles limits and vlues
            limperc = np.array([2.5, 97.5])    
            percen = np.linspace(limperc[0],limperc[1],num=self.nCells)
            
            #compute the 95% of the values from the population response
            for j in range(0,len(expTime)):
                populationlim[:,j] = np.percentile(AllcellExpr[:,j,idxPlot], percen[:])
            #end for j
            
            #compute minumum and maximum values in each instant of time
            ExprMin = np.amin(populationlim, axis=0)
            ExprMax = np.amax(populationlim, axis=0)
            
            #plot population output
            self.graphPlot3.canvas.axes.plot(expTime, np.median(AllcellExpr[:,:,idxPlot], axis=0), label=self.SpStr[idxPlot], linewidth=2)
            self.graphPlot3.canvas.axes.legend(loc='best')
            self.graphPlot3.canvas.axes.fill_between(expTime, ExprMin,ExprMax, alpha = 0.4)
            
            #highlight negative values in the plot
            if np.amin(AllcellExpr[:,:,idxPlot]) < 0:
                self.graphPlot3.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,idxPlot]),
                                                 color = 'C3', alpha = 0.1)
                self.graphPlot.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,idxPlot]),
                                                 color = 'C3', alpha = 0.1)
            
            #plot secondary population species
            if self.radioPlot2.isChecked() and self.subplotSP == True:
                #since the obtained values it takes the 95 % of the contained responses
                #between the 2.5 and 97.5 percentiles
                populationlim2 = np.zeros((AllcellExpr[:,:,self.indexPlot2].shape))
                
                #calculo del 95% de los valores de las respuestas de la población
                #compute the 95% values of the population output
                for j in range(0,len(expTime)):
                    populationlim2[:,j] = np.percentile(AllcellExpr[:,j,self.indexPlot2], percen[:])
                #end for j

                #compute minimum and maximum values in each instant of time
                ExprMin2 = np.amin(populationlim2, axis=0)
                ExprMax2 = np.amax(populationlim2, axis=0)   

                #plot population output 
                self.graphPlot5.canvas.axes.plot(expTime, np.median(AllcellExpr[:,:,self.indexPlot2], axis=0), label=self.SpStr[self.indexPlot2], linewidth=2)
                self.graphPlot5.canvas.axes.legend(loc='best')
                self.graphPlot5.canvas.axes.fill_between(expTime, ExprMin2,ExprMax2, alpha = 0.4)                 
                
                #highligh negative values in the plots
                if np.amin(AllcellExpr[:,:,:]) < 0:
                    self.graphPlot4.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,self.indexPlot2]),
                                                     color = 'C3', alpha = 0.1)
                    self.graphPlot5.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,self.indexPlot2]),
                                                     color = 'C3', alpha = 0.1)
        else:
            #plot when there is not variability
            self.graphPlot.canvas.axes.plot(expTime, AllcellExpr[idxPlot,:], label=self.SpStr[idxPlot])
            self.graphPlot.canvas.axes.legend(loc='best')
            
            #plot population output
            self.graphPlot3.canvas.axes.plot(expTime, AllcellExpr[idxPlot,:], linewidth=2, label=self.SpStr[idxPlot])
            self.graphPlot3.canvas.axes.legend(loc='best')
            
            #secondary species plot
            if self.radioPlot2.isChecked() and self.subplotSP == True:
                self.graphPlot4.canvas.axes.plot(expTime, AllcellExpr[self.indexPlot2,:], label=self.SpStr[self.indexPlot2])
                self.graphPlot4.canvas.axes.legend(loc='best')
                
                #plot population output
                self.graphPlot5.canvas.axes.plot(expTime, AllcellExpr[self.indexPlot2,:], linewidth=2, label=self.SpStr[self.indexPlot2])
                self.graphPlot5.canvas.axes.legend(loc='best')

        #draw plots into objects
        self.graphPlot.canvas.draw()
        self.graphPlot2.canvas.draw()
        self.graphPlot3.canvas.draw()
        self.graphPlot4.canvas.draw()        
        self.graphPlot5.canvas.draw()        
        
        #disable secondary species plot
        self.subplotSP = False
        
        self.labeliter.setText("Plot Done!")
        QApplication.processEvents()
  
    def computeSim(self, nCells, regressor): 
        #tag to show simulation state
        self.labeliter.setText("Simulating...")
        QApplication.processEvents()

        #array to save individual outputs
        AllcellExpr = np.zeros((nCells, len(regressor["inpU"]), len(regressor["ODEs"])))

        #system simulation depending on the variability type
        if self.noiseBox.isChecked() or self.intrinsicBox.isChecked() or self.extrinsicBox.isChecked():
            #process extrinsic variability
            if self.extrinsicBox.isChecked():
                print("MxEf")
                #compute natural log from parameter values
                self.logVpars = np.log(regressor["vPars"])

                #create log-covariance matrix
#                cov = 0.05*(np.identity(len(regressor["vPars"])))
                self.covStr = self.covMatrix.toPlainText()
                self.covStr = self.covStr.split(';')
                self.covStr = [x.split(',') for x in self.covStr]
                
                #convert logcovariance matrix
                self.valCov = list()
                for i in self.covStr:
                    self.valCov.append(list(map(float, i)))
                self.valCov = np.array(self.valCov)
                print(self.valCov)
                
                #check log-covariance matrix dimentions to be equal to number of parameters
                if self.valCov.shape[0] != len(regressor["vPars"]) or self.valCov.shape[1] != len(regressor["vPars"]):
                    self.outODEs.append("\nThe dimensions of the covariance matrix are not equal to the number of kinetic parameters.") 
                    QApplication.processEvents()
                    return
                
                #check log-covariance matrix to be symmetric
                elif not(self.isSymmetric(self.valCov, len(self.logVpars))):
                    self.outODEs.append("\nThe covariance matrix must be symmetric.") 
                    QApplication.processEvents()
                    return
                
                #compute cell multiparameters
                self.multipars = np.exp(np.random.multivariate_normal(self.logVpars, self.valCov, nCells))
                print(self.logVpars)
                print(self.valCov)
                print(self.multipars)

            else:
                #array with original parameters
                self.multipars = np.tile(regressor["vPars"], (nCells,1))
            print(self.multipars)
            
            #population simulations
            for k in range(0,nCells):
                print('Célula: ' + str(k+1))
                #tag to show simulation iteration
                self.labeliter.setText("Simulating... " + str(k+1))
                QApplication.processEvents()
            
                #process intrinsic variability
                if self.intrinsicBox.isChecked():
                    print("Gill")
                    #numeric solution of ODEs
                    try:
                        cellExpr = s2b.cGill(self.multipars[k,:], regressor)
                    except:
                        cellExpr = s2b.gillAl(self.multipars[k,:], regressor)
                else:
                    #numeric solution of ODEs
                    cellExpr = s2b.solveODE(self.multipars[k,:], regressor)
                    
                #store individual simulated outputs
                for j in range(0, len(cellExpr)):
                    AllcellExpr[k,:,j] = cellExpr[j,:]
            
            if self.noiseBox.isChecked():
                print("Noise")
                
                #get noise parameters
                self.changeErrorPar()
                
                #normal distribution
                mu, sigma = 0, 1 #mean and standard deviation
                ndis = np.random.normal(mu, sigma, (nCells,len(regressor["inpU"])))
                
                #compute individual output with noise for each species
                hvar = self.a_err + self.b_err*AllcellExpr[:,:,-1]     
                noisySp = AllcellExpr[:,:,-1] + hvar*ndis
                
                #store species profiles with noise
                AllcellExpr[:,:,-1] = noisySp
                
        else:
            #deterministic output
            #numeric solution of ODEs
            cellExpr = s2b.solveODE(regressor["vPars"], regressor)
            AllcellExpr = cellExpr
        
        print(self.noiseBox.isChecked(), self.intrinsicBox.isChecked(), self.extrinsicBox.isChecked())
                
        self.labeliter.setText("Simulation Done!")
        QApplication.processEvents()
        
        return AllcellExpr
    
    #function to verify if log-covariance matrix is symmetric
    def isSymmetric(self, mat, N): 
        for i in range(N): 
            for j in range(N): 
                #compare matrix values with its equivalent
                if (mat[i,j] != mat[j,i]): 
                    return False
        return True
    
    #choose input type
    def changeInput(self):
        self.inputT = self.inputType.currentText()
        
        if self.inputT == 'Model':
            self.setButton.setEnabled(True)
        else:
            self.setButton.setEnabled(False)
        
    #show simulation information
    def showInfosim(self, expTdur, concenStr, nCells):
        #show experiment duration
        self.outODEs.append("\nDuration of the Experiment: " + str(expTdur) + " min")
        
        #show initial concentration values
        SpConcent = [(self.SpStr[i] + ' = ' + concenStr[i]) for i in range(len(concenStr))]
        print(SpConcent)
        self.outODEs.append("Initial Concentration: " + ', '.join(SpConcent))  
        
        if self.idxReaction != None:
            pulseStr = [(self.tOnStr[i] + " (" + self.tDurStr[i] + ")") for i in range(len(self.valTon))]
            #show input pulses
            self.outODEs.append("Input Pulses: " + ', '.join(pulseStr))
            
        #show number of cells
        self.outODEs.append("Number of Cells: " + str(nCells))

    #function to update experiment duration value
    def changeTdur(self):
        self.expTdur = self.spinTdur.value()
        
    #get noise parameters
    def changeErrorPar(self):
        self.a_err = self.spinAval.value()
        self.b_err = self.spinBval.value()
        print(self.a_err, self.b_err)
        
    #function to update number of cells
    def changeCells(self):
        self.nCells = self.numCells.value()
    
    #open window to define input model
    def setInputModel(self):
        self.inputModelW.exec_()
    
    #enable input box to define log-covarince matrix
    def enable_cov(self):
        if self.extrinsicBox.isChecked():
            self.covMatrix.setEnabled(True)
        else:
            self.covMatrix.setEnabled(False)
    
    #enable input to define noise parameter values
    def enable_noise(self):
        if self.noiseBox.isChecked():
            self.spinAval.setEnabled(True)
            self.spinBval.setEnabled(True)
        else:
            self.spinAval.setEnabled(False)
            self.spinBval.setEnabled(False)
               
    #change species to be plotted
    def changePlot(self):
        #clean plots
        self.clearPlot()
        #disable while plitting
        self.changePlotButton.setEnabled(False)
                
        #get species to be plotted
        print(self.spPlotBox.currentText())
        self.idxPlot2 = self.SpStr2.index(self.spPlotBox.currentText())
        print(self.idxPlot2)
        
        #get secondary species to be plotted
        if self.radioPlot2.isChecked():
            self.indexPlot2 = self.SpStr.index(self.spPlotBox2.currentText())
            print(self.indexPlot2)
            print(self.spPlotBox2.currentText())
            #enable secondary species subplot
            self.subplotSP = True
        
        #plot chose species
        if not(self.spPlotBox.currentText() == "All"):
            self.plotting(self.expTime, self.inpShock, self.AllcellExpr, self.idxPlot2)
        
        else:
            #plot all outputs
            for i in range(len(self.SpStr)):
                self.plotting(self.expTime, self.inpShock, self.AllcellExpr, i)
        
        #enable after plotting
        self.changePlotButton.setEnabled(True)
        #disable secondary species subplot
        self.subplotSP = False
    
    #enable combobox to choose secondary species to be plotted
    def enablePlot2(self):
        if self.radioPlot2.isChecked():
            self.spPlotBox2.setEnabled(True)
        else:
            self.spPlotBox2.setEnabled(False)
                    
    def inferate(self):
        #clean plots
        self.clearPlot()
        #restore size of main plot
        self.graphPlot.resize(891,411)

        #disable while infering
        self.enableSimOpt(False, False)

        #set tab 1
        self.tabWidget.setCurrentIndex(0)

        #get type modelo to perform inference
        self.modelInfr = self.modelBox.currentText()
        print(self.modelInfr)
        
        #take initial parameter values
        self.iparStr = self.pars0.text()
        self.iparStr = re.findall('[0-9.0-9]+', self.iparStr)
        
        try:
            self.beta0 = list(map(float, self.iparStr))
            if len(self.beta0) != len(self.regressor["vPars"]):
                self.outODEs.append("The number of initial parameters must be equal to the number of kinetic parameters.")
                self.enableSimOpt(True, True)
                return
            print(self.beta0)
            
            if len(self.data["genExpr"].shape) == 3:
                self.Obs = self.data["genExpr"][:,:,-1]
                print(self.Obs.shape)
            elif len(self.data["genExpr"].shape) == 2:
                self.Obs = np.array([self.data["genExpr"][-1,:]])
                print(self.Obs.shape)
            
            self.betacal = self.computeInfr(self.modelInfr, self.beta0, self.Obs, self.regressor)
            
            self.showInfoInfr(self.betacal)
            
            #enable when inference is done
            self.enableSimOpt(True, True)
        
        except:
            self.enableSimOpt(True, True)
            self.outODEs.append("\nPlease review the input information.")

    #perform model stimations
    def computeInfr(self, model, beta0, Obs, regressor):
        self.labeliter.setText("Inferring...")
        QApplication.processEvents()                   
        
        #determine mean and standar deviation curves
        self.uObs = np.mean(Obs, axis=0)
        self.sdObs = np.std(Obs, axis=0)
        
        #initial noise parameter values
        self.A0 = self.A0Box.value()
        self.B0 = self.B0Box.value()
        self.errAB = np.array([self.A0, self.B0])
        print(self.errAB)
        
        #number of iterations to be perform
        self.iterInfr = self.iterBox.value()
        
        self.step = 0
        if model == "Mean Curve Fit":
            #use curve_fit functon to perform regression over the mean curve
            #to infer parameters
            def ODEfit(regressor, *pars):
                cellExpr = s2b.solveODE(pars, regressor)
                
                self.step += 1
                if self.step % 2:
                    print(pars)
                    self.plotInfr(self.expTime, self.uObs, cellExpr[-1,:])
                return cellExpr[-1,:]
            
            betacal, cov = curve_fit(ODEfit, regressor, self.uObs, beta0)
            print(betacal, len(betacal))
        
        elif model == "Mean Cell": 
            #use likelihood function to estimate population parameters
            regressor["allObs"] = Obs
            regressor["errPars"] = self.errAB

            self.step = 0
            def MLL_mc(pars, regressor):
                #compute simulation with temporal parameters
                cellExpr= s2b.solveODE(pars, regressor)
                meanP = cellExpr[-1,:]
                
                sdP = regressor["errPars"][0] + regressor["errPars"][1]*meanP
                sdP = np.absolute(sdP)
                
                self.step += 1
                if self.step % 2:
                    print(pars)
                    self.plotInfr(self.expTime, self.uObs, meanP, self.sdObs, sdP)
            
                ncells = len(regressor["allObs"])
                
                fSim = np.tile(meanP, (ncells,1))
                hSim = np.tile(sdP, (ncells,1))
                
                Observaciones = regressor["allObs"]
                
                mllmc = s2b.MLLmeasure(Observaciones, fSim, hSim)
                
                return mllmc
            
            #estimation of kinetics parameters
            betacal = fmin(MLL_mc, x0 = beta0, args=(regressor,),
                           maxiter=self.iterInfr, maxfun=self.iterInfr)
            print(betacal)
            regressor["betaCal"] = betacal
            
            #mean curve simulation from infered kinetic parameters
            cellBeta = s2b.solveODE(betacal, regressor)
            regressor["uC"] = cellBeta[-1,:]
            
            #minus log-likelihood function
            self.step = 0
            def MLL_mcErr(error, regressor):
                meanP = regressor["uC"]
                sdP = error[0] + error[1]*meanP
                sdP = np.absolute(sdP)
                
                self.step += 1
                if self.step % 2:
                    print(error)
                    self.plotInfr(self.expTime, self.uObs, meanP, self.sdObs, sdP)

                #change of cost value to fit model
                ncells = len(regressor["allObs"])
                fSim = np.tile(meanP, (ncells,1))
                hSim = np.tile(sdP, (ncells,1))
                
                #observations array
                Observaciones = regressor["allObs"]
                
                #compute minus-log-likelihood
                mllmc = s2b.MLLmeasure(Observaciones, fSim, hSim)
                return mllmc
            
            #minimization of cost function based on likelihood estimaton of
            #noise parameters
            result = fmin(MLL_mcErr, x0 = self.errAB, args=(regressor,),
                          maxiter=self.iterInfr, maxfun=self.iterInfr)
            self.noiseCal = result
            print(result)
            
            #plot curve with infered parameters
            uInfr= s2b.solveODE(betacal, regressor)
            sdInfr = self.noiseCal[0] + self.noiseCal[1]*uInfr[-1,:]
            self.plotInfr(self.expTime, self.uObs, uInfr[-1,:], self.sdObs, sdInfr)

        elif model == "Moment Based":
            #use KLD to fit moments and estimate parameters
          
            #compute moments
            #determine system without input
            if self.idxReaction == None:
                self.ODE2M, self.vars2M = s2b.simbMoments(self.SpStr,
                    self.valReac, self.valProd,self.nparStr)
                
            #determine system with input
            else:
                self.ODE2M, self.vars2M = s2b.simbMoments(self.SpStr,
                    self.valReac, self.valProd, self.nparStr, 
                    inputN=self.nameInp, indexU = self.idxReaction)
            print(self.ODE2M)
            
            #moment initial conditions
#            sp0 = np.concatenate((regressor["species0"], np.zeros(len(self.ODE2M) - len(regressor["species"]))))
            MomentSp0 = sym.lambdify([self.SpStr], self.vars2M['species'], "numpy")
            sp0M = np.array(MomentSp0(regressor["species0"]))
            print('Concentración inicial momentos', sp0M)
            
            #regressor to compute moments
            regressor2 = {
                    "ODEs": self.ODE2M,
                    "regressor":regressor,
                    "inpU": self.inpShock,
                    "Vtime": self.expTime,
                    "species0":sp0M,
                    "meanCell":self.uObs,
                    "sdCell":self.sdObs
                    }
            regressor2.update(self.vars2M)
            
            #find second moment index of target species
#            idx2M = list(map(str, regressor2["species"]))
#            idx2M = idx2M.index(self.SpStr[-1] + '**2')
#            print(idx2M)
#            regressor2["idx2M"] = idx2M
            regressor2["idx2M"] = -1
            regressor2["errPars"] = self.errAB
            
            #KLD cost function  
            self.step = 0
            def KLDmb(pars, data):
            
                uC, sdC = s2b.solve2Moment(pars, data["errPars"], data)
                
                self.step += 1
                if self.step % 2:
                    print(pars)
                    self.plotInfr(self.expTime, self.uObs, uC, self.sdObs, sdC)
                
                uObs = data["meanCell"]
                sdObs = data["sdCell"]
                
                kld = s2b.kldmeasure(uObs, sdObs, uC, sdC)
                ukld = np.mean(kld)
                return ukld
            
            #kinetic parameters estimations
            betacal = fmin(KLDmb, x0 = beta0, args=(regressor2,), 
                           maxiter=self.iterInfr, maxfun=self.iterInfr)
            print(betacal)
            regressor["betaCal"] = betacal
            
            #compute moments
            momentCell = s2b.solveODE(betacal, regressor2)
            regressor2["uC"] = momentCell[len(regressor["species"]) - 1,:]
            regressor2["sd2C"] = momentCell[regressor2["idx2M"],:]
            
            #cost function to minimize according to kld measure
            self.step = 0
            def KLDmbErr(errorP, data):
                print(errorP)
                uC = data["uC"]
                sd2C = data["sd2C"] - uC**2
                
                #compute standar deviation of the second order moment from
                #raw variance 
                sdC = np.sqrt((1 + errorP[1]**2)*sd2C + (errorP[0] + errorP[1]*uC)**2)
                
                #observed data
                uObs = data["meanCell"]
                sdObs = data["sdCell"]
                
                self.step += 1
                if self.step % 2:
                    print(errorP)
                    self.plotInfr(self.expTime, self.uObs, uC, self.sdObs, sdC)
                
                kld = s2b.kldmeasure(uObs, sdObs, uC, sdC)
                ukld = np.mean(kld)
                return ukld
            
            #compute parametes through kld
            result = fmin(KLDmbErr, x0 = self.errAB, args=(regressor2,), 
                          maxiter=self.iterInfr, maxfun=self.iterInfr)
            self.noiseCal = result
            print(self.noiseCal)
            
            #plot curve with infered parameters
            uC, sdC = s2b.solve2Moment(betacal, self.noiseCal, regressor2)
            self.plotInfr(self.expTime, self.uObs, uC, self.sdObs, sdC)
        
        self.labeliter.setText("Inference Done!")
        
        return betacal
    
    #plot infered curve in each iteration
    def plotInfr(self, expTime, uCell, InfrCell, sdCell=np.array([]), sdInfrCell=np.array([])):
        self.graphPlot.canvas.axes.clear()
        self.graphPlot.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot.canvas.axes.plot(expTime, uCell, label="Observed Output")
        self.graphPlot.canvas.axes.plot(expTime, InfrCell, label="Inferred Output")
        self.graphPlot.canvas.axes.set_ylim(np.amin(uCell), np.amax(uCell))
        self.graphPlot.canvas.axes.legend(loc='best')
        self.graphPlot.canvas.axes.grid()
        
        #plotting of output deviation
        if sdCell.size != 0:
            downL = uCell - sdCell
            upL = uCell + sdCell
            self.graphPlot.canvas.axes.plot(expTime, downL, 'b-')
            self.graphPlot.canvas.axes.plot(expTime, upL, 'b-')
            self.graphPlot.canvas.axes.plot(expTime, InfrCell - sdInfrCell, 'y-')
            self.graphPlot.canvas.axes.plot(expTime, InfrCell + sdInfrCell, 'y-')
            self.graphPlot.canvas.axes.set_ylim(np.amin(downL), np.amax(upL))
        
        self.graphPlot.canvas.draw()
        QApplication.processEvents()                   
    
    #show infered information
    def showInfoInfr(self, betacal):
        
        #show moments
        if self.modelInfr == "Moment Based":
            self.outODEs.append("\nDifferential Equations of Moments: ")
            for j in range(len(self.vars2M["spODE"])):
                self.outODEs.append('\t' + 'd' + str(self.vars2M["species"][j]) + '/dt =  ' + str(self.vars2M["spODE"][j]))
        
        #show infered parameters
        nparStrV = [(self.nparStr[i] + ' = ' + str(betacal[i])) for i in range(len(betacal))]
        self.outODEs.append("\nInferred Parameters:\n " + ', '.join(nparStrV))
        
        #show noise parameters
        if self.modelInfr != "Mean Curve Fit":
            self.outODEs.append("\na = " + str(self.noiseCal[0]) + ', b = ' + str(self.noiseCal[1]))
    
    #function to change inference model
    def changeModel(self):
        if self.modelBox.currentText() == "Mean Curve Fit":
            self.A0Box.setEnabled(False)
            self.B0Box.setEnabled(False)
            self.iterBox.setEnabled(False)
        else:
            self.A0Box.setEnabled(True)
            self.B0Box.setEnabled(True)
            self.iterBox.setEnabled(True)

    #start and end events of window execution
    def showEvent(self, event):
        self.outODEs.setText("Welcome to GenExpSim!")
    
    def closeEvent(self, event):
        pushClose = QMessageBox.question(self, "Exit",
                                         "Are you sure you want to exit??",
                                         QMessageBox.Yes | QMessageBox.No)
        if pushClose == QMessageBox.Yes: event.accept()
        else: event.ignore() 
        
#instance to start application
app = QApplication(sys.argv)
#object of the created class
_mainWindow = mainWindow()
#show window
_mainWindow.show()
#execute application
app.exec_()