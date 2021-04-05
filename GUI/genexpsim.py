# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 07:45:55 2020

@author: Bry
"""
#imports libraries used to programm GUI
import sys, re
import os
import traceback
from datetime import datetime
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QDesktopWidget, QDialog, QAction, QFileDialog, QInputDialog
from PyQt5.QtGui import QIcon
#package to load the .ui file
from PyQt5 import uic
#libraries for processing data
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import fmin
import cma
import sympy as sym
import simsysbio as s2b
import time 

#builds the dialog window used to define the input model
class inputWindow(QDialog):
    def __init__(self):
        QDialog.__init__(self)
        uic.loadUi("inputWindow.ui", self)
        
        #connects button  to define input model
        self.defineInputButton.clicked.connect(self.defineInput)
        #validates information from inputs
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
        
        try:
            #takes parameter values
            parStrInp = self.inpModelPar.text()
            parStrInp = re.findall('[0-9.0-9]+', parStrInp)
            valParInp = list(map(float, parStrInp))
            self.inpReacBox.setRange(1, len(valParInp))
        except:
            pass
    
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
    
    #takes data to define input model when "Define" is clicked
    def defineInput(self):
        #takes the name of the molecular species
        self.SpStrInp = self.inpModelSp.text()
        self.SpStrInp = self.SpStrInp.split(',')
        self.SpStrInp = [x.strip() for x in self.SpStrInp]
        
        #takes parameter values
        self.parStrInp = self.inpModelPar.text()
        self.parStrInp = re.findall('[0-9.0-9]+', self.parStrInp)
        
        #takes reagent matrix
        self.reacStrInp = self.inpModelR.toPlainText()
        self.reacStrInp = self.reacStrInp.split(';')
        self.reacStrInp = [x.split(',') for x in self.reacStrInp]
        
        #takes product matrix
        self.prodStrInp = self.inpModelP.toPlainText()
        self.prodStrInp = self.prodStrInp.split(';')
        self.prodStrInp = [x.split(',') for x in self.prodStrInp]
        
        try:
            #converts parameters and get the name from each one
            self.valParInp = list(map(float, self.parStrInp))
            self.nparStrInp = [('k'+ str(i+1)) for i in range(len(self.parStrInp))]
            
            #converts data from reagent matrix
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
            
            self.inpindex = self.inpReacBox.value()
            
            #compute the simbolic expression of the input model
            self.ecuacionesInp, self.variablesInp = s2b.simbODE(self.SpStrInp,
                        self.valReacInp, self.valProdInp, self.nparStrInp, 
                        inputN=self.nameInp, indexU=self.inpindex)
            
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
            message = traceback.format_exc()
            self.msgBox(message)    

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
        #create window instance of the model input definition
        self.inputModelW = inputWindow()
        #create toolbar
        self.Create_toolbar()
        
        #initial state of some variables
        self.idxReaction = None
        #dictionary to save data that plays as regressor
        self.regressor = dict()
        #disable inputs that depend on other options
        self.setButton.setEnabled(False)
        self.spinAval.setEnabled(False)
        self.spinBval.setEnabled(False)
        self.covMatrix.setEnabled(False)
        self.changePlotButton.setEnabled(False)
        self.sampleCells.setEnabled(False)
        #initial state of the secondary species subplot
        self.subplotSP = False
        self.a_err = 0.0
        self.b_err = 0.0
        self.tempAxes = []
        self.valTon = [0]
        self.tOnStr = [0]
        self.valTdur = [0]
        self.tDurStr = [0]
                
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
        
    def Create_toolbar(self):
        saveAct = QAction(QIcon('save.png'), 'Save', self)
        saveAct.setShortcut('Ctrl+S')
        saveAct.triggered.connect(self.saveData)
        
        openAct = QAction(QIcon('open.png'), 'Open', self)
        openAct.setShortcut('Ctrl+O')
        openAct.triggered.connect(self.openData)
        
        self.toolbar = self.addToolBar('Toolbar')
        self.toolbar.addAction(saveAct)
        self.toolbar.addAction(openAct)
        
    def saveData(self):
        #show status
        self.labeliter.setText("Saving Data...")
        QApplication.processEvents()
        
        #get the path to the folder where data is stored
        nowPath = os.getcwd() + '/data'
        #create a folder to save data
        try:
            os.makedirs('data/')
        except:
            pass
        
        try:
            #open window explorer to save data
            name = QFileDialog.getSaveFileName(self, 'Save Data', nowPath, 
                                               "Numpy File (*.npz)")
            
            #save data into a .npz file
            np.savez_compressed(name[0], SpStr=self.SpStr,
                     parStr=self.parStr,
                     valReac=self.valReac,
                     valProd=self.valProd,
                     ecuaciones=list(map(str, self.ecuaciones)),
                     vPars=np.array(self.valPar),
                     idxR=self.idxReaction,
                     species=list(map(str, self.variables['species'])),
                     pars=list(map(str, self.variables['pars'])),
                     nameVar=list(map(str, self.variables['nameVar'])),
                     species0=np.array(self.valConcen),
                     Vtime=self.expTime,
                     inpU=self.inpShock,
                     nCells=self.nCells,
                     concenStr=self.concenStr,
                     expTdur=self.expTdur,
                     valTon=self.valTon,
                     tOnStr=self.tOnStr,
                     valTdur=self.valTdur,
                     tDurStr=self.tDurStr,
                     AllcellExpr=self.AllcellExpr,
                     nparStr=self.nparStr,
                     nameInp=self.nameInp,
                     a_err=self.a_err,
                     b_err=self.b_err)
            
            self.labeliter.setText("Data Saved!")
            
        except:
            self.labeliter.setText("")
            if not name[0] == '':
                message = traceback.format_exc()
                self.msgBox(message)
        
        print(name[0])
        print('Data saved')
        
    def openData(self):
        self.labeliter.setText("Loading Data...")
        QApplication.processEvents()
        
        #get the path to the folder where data is stored
        nowPath = os.getcwd() + '/data'
        
        try:
            #open window explorer to find stored data
            name = QFileDialog.getOpenFileName(self, 'Open Data', nowPath, 
                                               "Numpy File (*.npz)")
            
            #load data from .npz file
            npzfile = np.load(name[0])
            
            #get data saved
            #shows data added for defining
            self.SpStr = npzfile['SpStr'].tolist()
            self.parStr = npzfile['parStr'].tolist()
            self.valReac = npzfile['valReac']
            self.valProd = npzfile['valProd']
            self.ecuaciones = npzfile['ecuaciones'].tolist()
            self.valPar = npzfile['vPars']
            self.expTime = npzfile['Vtime']
            self.inpShock = npzfile['inpU']
            self.nCells = npzfile['nCells']
            self.concenStr = npzfile['concenStr']
            self.expTdur = npzfile['expTdur']
            self.valTon = npzfile['valTon']
            self.tOnStr = npzfile['tOnStr']
            self.valTdur = npzfile['valTdur']
            self.tDurStr = npzfile['tDurStr']
            self.AllcellExpr = npzfile['AllcellExpr']
            self.nparStr = npzfile['nparStr'].tolist()
            self.nameInp = str(npzfile['nameInp'])
            self.valConcen = npzfile['species0'].tolist()
            self.a_err = npzfile['a_err']
            self.b_err = npzfile['b_err']
            
            try:
                self.idxReaction = npzfile['idxR']
            except:
                self.idxReaction = None
                
            #computes system without input
            if self.idxReaction == None:
                self.ecuaciones, self.variables = s2b.simbODE(self.SpStr,
                    self.valReac, self.valProd,self.nparStr)
                
            #Computes system with input
            else:
                self.ecuaciones, self.variables = s2b.simbODE(self.SpStr,
                    self.valReac, self.valProd, self.nparStr, 
                    inputN=self.nameInp, indexU = self.idxReaction)
            
            #stores information used in simulation
            self.regressor = {
                    "ODEs": self.ecuaciones,
                    "matrizR": self.valReac,
                    "matrizP": self.valProd,
                    "vPars": self.valPar,
                    "idxR": self.idxReaction,
                    "species0":npzfile['species0'],
                    "Vtime":npzfile['Vtime'],
                    "inpU":npzfile['inpU']
                    }
            
            self.variables = {
                    "species":npzfile['species'],
                    "pars":npzfile['pars'],
                    "nameVar":npzfile['nameVar']
                    }
            self.regressor.update(self.variables)
            
            ##### Identify loaded section ####
            self.showInfo(self.SpStr, self.parStr, self.valReac, self.valProd, 
                          self.ecuaciones)
            
            #creates .py file containig the differential equations system
            s2b.modelDefiner(self.variables["species"], self.ecuaciones, 
                             self.variables["pars"])
            
            #add items to change species plot
            self.addItemsPlot()
            
            #### #simulate loaded section ####
            self.sampleCells.setMaximum(self.nCells)
            
            #shows information for simulating
            self.showInfosim(self.expTdur, self.concenStr, self.nCells)
            
            #stores simulated data
            self.data = {
                    "genExpr": self.AllcellExpr,
                    "inPulse": self.inpShock,
                    "vTime": self.expTime,
                    }
            
            #plots system output
            self.clearPlot()
            self.idxPlot = self.SpStr.index(self.outPlot.currentText())
            self.plotting(self.expTime, self.inpShock, self.AllcellExpr, 
                          self.idxPlot)
            
            #set plot button and add items
            self.setPlotButton()
    
            #actives button: Simulation
            self.enableSimOpt(True, True)
            
            self.labeliter.setText("Data Loaded!")
            print(npzfile.files)
            
        except:
            self.labeliter.setText("")
            if not name[0] == '':
                message = traceback.format_exc()
                self.msgBox(message)
        
        print(name[0])
        print('Data opened')
    
    def plotSetting(self):
        self.graphPlot.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot.canvas.axes.grid()
        
        self.graphPlot4.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot4.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot4.canvas.axes.grid()
                
        self.graphPlot2.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot2.canvas.axes.set_ylim([0, 1])
        self.graphPlot2.canvas.axes.set_xlabel("Time (min)", labelpad=0.0, 
                                               fontsize=8)
        self.graphPlot2.canvas.axes.grid()
        
        self.graphPlot3.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot3.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot3.canvas.axes.grid()

        self.graphPlot5.canvas.axes.set_ylabel("Concentration", fontsize=8)
        self.graphPlot5.canvas.axes.tick_params(bottom=False, labelbottom=False)
        self.graphPlot5.canvas.axes.grid()
        
        self.graphPlot6.canvas.axes.set_title("Population Parameter Iteration",
                                              fontsize=10)        
        self.graphPlot6.canvas.axes.tick_params(bottom=False, left=False, 
                                                labelbottom=False, 
                                                labelleft=False)
        self.graphPlot6.canvas.axes.set_xlabel("Iterations", labelpad=16.0, 
                                               fontsize=8)
        
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
        #enables inputs of variability sources
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
        #checks input information
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
        #combobox to choose optimization algorithm
        self.AlgoBox.currentIndexChanged.connect(self.changeAlgorithm)
        #change to parameters tab
        self.tabWidget.currentChanged.connect(self.tabChange)
        
    #change tab size when parameters chart tab is selected.
    def tabChange(self):
        if self.tabWidget.currentIndex() == 2:
            self.tabWidget.resize(781,528)
        else:
            self.tabWidget.resize(781,371)
        
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
    
    #clears all plots
    def clearPlot(self):
        self.graphPlot.canvas.axes.clear()
        self.graphPlot2.canvas.axes.clear()
        self.graphPlot3.canvas.axes.clear()
        self.graphPlot4.canvas.axes.clear()
        self.graphPlot5.canvas.axes.clear()
        self.graphPlot6.canvas.axes.clear()
        
        self.clearSubplot()
            
        self.plotSetting()

        self.graphPlot.canvas.draw()
        self.graphPlot2.canvas.draw()
        self.graphPlot3.canvas.draw()
        self.graphPlot4.canvas.draw()
        self.graphPlot5.canvas.draw()
        self.graphPlot6.canvas.draw()
        
        self.tabWidget.setTabText(0, "Individual")
        self.tabWidget.setTabText(1, "Population")
    
    def clearSubplot(self):
        if not len(self.tempAxes) == 0:
            for i in range(len(self.tempAxes)):
                self.tempAxes[i].remove()
            self.tempAxes = []

    #enables simulation and inference section
    def enableSimOpt(self, eSim = False, eOpt = False):
        #enables simulation
        self.simButton.setEnabled(eSim)
        self.tabSim.setEnabled(eSim)
        
        #Enables input configuration
        if self.idxReaction == None:
            self.inputType.setEnabled(False)
            self.tOnPulse.setEnabled(False)
            self.tDurPulse.setEnabled(False)
        else:
            self.inputType.setEnabled(True)
            self.tOnPulse.setEnabled(True)
            self.tDurPulse.setEnabled(True)
        
        #enables inference
        self.optButton.setEnabled(eOpt)
        self.tabOpt.setEnabled(eOpt)
        
    #action executed when "define" is clicked
    def identify(self):
        #disables other buttons or clean inputs
        self.enableSimOpt()
        self.labeliter.clear()
        self.clearPlot()
        
        #shows definition state
        self.labeliter.setText("Definig System...")
        QApplication.processEvents()
        
        #takes molecular species names
        self.SpStr = self.inEspecies.text()
        self.SpStr = self.SpStr.split(',')
        self.SpStr = [x.strip() for x in self.SpStr]
        
        #takes parameter values
        self.parStr = self.inPars.text()
        self.parStr = re.findall('[0-9.0-9]+', self.parStr)
        
        #takes reagent matix values
        self.reacStr = self.inReact.toPlainText()
        self.reacStr = self.reacStr.split(';')
        self.reacStr = [x.split(',') for x in self.reacStr]
        
        #takes product matrix values
        self.prodStr = self.inProd.toPlainText()
        self.prodStr = self.prodStr.split(';')
        self.prodStr = [x.split(',') for x in self.prodStr]
                
        try:
            #shows data added for definig system
            self.showInfo(self.SpStr, self.parStr)
            
            #converts parameters and get parameter names
            self.valPar = list(map(float, self.parStr))
            self.nparStr = [('c'+ str(i+1)) for i in range(len(self.parStr))]
            print(self.valPar)
            print(self.nparStr)
            print(self.SpStr)
            
            #converts data from reagent matrix
            self.valReac = list()
            for i in self.reacStr:
                self.valReac.append(list(map(int, i)))
            self.valReac = np.array(self.valReac)
            print(self.valReac)
            
            #converts data from product matrix
            self.valProd = list()
            for i in self.prodStr:
                self.valProd.append(list(map(int, i)))
            self.valProd = np.array(self.valProd)
            print(self.valProd)
            print(self.idxReaction)
            
            #checks input dimentions
            if self.valReac.shape != self.valProd.shape:
                self.outODEs.append("\nThe dimensions of the matrices are not equal.") 
                return
            elif (len(self.valReac) != len(self.SpStr)) or (len(self.valProd) != len(self.SpStr)):
                self.outODEs.append("\nThe number of rows in the matrices is not equal to the number of species.") 
                return
            elif (self.valReac.shape[1] != len(self.valPar)) or (self.valProd.shape[1] != len(self.valPar)):
                self.outODEs.append("\nThe number of columns in the matrices is not equal to the number of parameters.") 
                return
            
            #computes simbolic expressions of the system
            #gets ODEs in simbolic form together to their respective varible names
            self.nameInp = 'u'
            
            #computes system without input
            if self.idxReaction == None:
                self.ecuaciones, self.variables = s2b.simbODE(self.SpStr,
                    self.valReac, self.valProd,self.nparStr)
                
            #Computes system with input
            else:
                self.ecuaciones, self.variables = s2b.simbODE(self.SpStr,
                    self.valReac, self.valProd, self.nparStr, 
                    inputN=self.nameInp, indexU = self.idxReaction)
            
            #creates .py file containig the differential equations system
            s2b.modelDefiner(self.variables["species"], self.ecuaciones, 
                             self.variables["pars"])

            #shows status
            self.labeliter.setText("System Defined!")

            #shows data added for defining
            self.showInfo(self.SpStr, self.parStr, self.valReac, self.valProd,
                          self.ecuaciones)
            
            #actives button: Simulation
            self.enableSimOpt(True)
            
            #stores information used in simulation
            self.regressor = {
                    "ODEs": self.ecuaciones,
                    "matrizR": self.valReac,
                    "matrizP": self.valProd,
                    "vPars": np.array(self.valPar),
                    "idxR": self.idxReaction
                    }
            self.regressor.update(self.variables)
            
            self.addItemsPlot()

            #noise measurement messagge
            self.noiseBox.setToolTip("Measurement noise only will affect \n the observed species: " + str(self.SpStr[-1]))
            
        except:
            self.outODEs.append("\nPlease review the input information.")
            message = traceback.format_exc()
            self.msgBox(message)
            
    def addItemsPlot(self):
        #set species names in combobox
        self.outPlot.clear()
        self.SpStr.reverse()
        self.outPlot.addItems(self.SpStr)
        self.SpStr.reverse()
            
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

            #shows differential equations
            self.outODEs.append("Differential Equations System: ")
            for j in range(len(ecuaciones)):
                self.outODEs.append('\t' + 'd' + SpStr[j] + '/dt =  ' + str(ecuaciones[j]))
                    
    #takes reaction-input index if it changes
    def changeIdxR(self):
        self.idxReaction = self.Ridx.value()

    def enable_entrada(self):
        #enables or disables inputs for definig of one system input
        self.parStr = None
        self.parStr = self.inPars.text()
        self.parStr = self.parStr.split(',')
        
        #enables input just if it detects any parameter
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
        #sets tab 1
        self.tabWidget.setCurrentIndex(0)
        #reboots current and previous sections
        self.clearPlot()
        self.changePlotButton.setEnabled(False)

        #disable during simulation
        self.enableSimOpt(False)
        
        #takes experiment duration value
        self.expTdur = self.spinTdur.value()
        print(self.expTdur)
        
        #number of cells to simulate
        self.nCells = self.numCells.value()
        self.sampleCells.setMaximum(self.nCells)
        
        #takes initial concentrations values
        self.concenStr = self.concentSp.text()
        self.concenStr = re.findall('[0-9.0-9]+', self.concenStr)
        
        #takes input configuration if there is one
        if self.idxReaction != None:
            #takes pulse start values
            self.tOnStr = self.tOnPulse.text()
            self.tOnStr = re.findall('[0-9.0-9]+', self.tOnStr)
        
            #takes pulse duration values
            self.tDurStr = self.tDurPulse.text()
            self.tDurStr = re.findall('[0-9.0-9]+', self.tDurStr)
        
        try:
            #converts initial concentration values to float
            self.valConcen = list(map(float, self.concenStr))
            self.regressor["species0"] = np.array(self.valConcen)
        
            #checks concentration values to be equal to amount of species
            if len(self.SpStr) != len(self.valConcen):
                self.outODEs.append("\nThe number of initial concentration values is not equal to the number of species.")
                self.enableSimOpt(True)
                return
            print(self.valConcen)
            
            #computes input signal
            if self.idxReaction == None:
                #time and input vectors
                self.expTime = np.linspace(0, self.expTdur, self.expTdur + 1)
                self.inpShock = np.ones(len(self.expTime))
                
            else:
                #convets pulse values to float
                self.valTon = list(map(float, self.tOnStr))
                self.valTdur = list(map(float, self.tDurStr))
                print(self.valTon, self.valTdur)
                
                #checks experiment duration
                if not(all((i >= 0 and i < self.expTdur) for i in self.valTon)):
                    self.outODEs.append("\nThe start of the pulses is not within the established limits.")
                    self.enableSimOpt(True)
                    return
                elif not(all((i > 0) for i in self.valTdur)):
                    self.outODEs.append("\nThe duration of the pulses must be greater than zero.")
                    self.enableSimOpt(True)
                    return
                
                #gets input type
                self.changeInput()
                print(self.inputT, type(self.inputT))
                
                #computes input system
                self.expTime, self.inpShock = self.computeInput()
            
            #saves information for simulating
            self.regressor["Vtime"] = self.expTime
            self.regressor["inpU"] = self.inpShock
            print(self.regressor)
            
            #shows information for simulating
            self.showInfosim(self.expTdur, self.concenStr, self.nCells)
            
            #computes system simulation
            self.AllcellExpr = self.computeSim(self.nCells, self.regressor)
            try:
                print(self.AllcellExpr.shape)
            except:
                print(self.AllcellExpr)
                
            #plots system output
            self.idxPlot = self.SpStr.index(self.outPlot.currentText())
            self.plotting(self.expTime, self.inpShock, self.AllcellExpr, self.idxPlot)
                        
            #stores simulated data
            self.data = {
                    "genExpr": self.AllcellExpr,
                    "inPulse": self.inpShock,
                    "vTime": self.expTime,
                    }
            
            #enables plot button and add items
            self.setPlotButton()
            
            #enables inference section
            self.enableSimOpt(True, True)
        except:
            self.enableSimOpt(True)
            self.outODEs.append("\nPlese, review the input information")
            message = traceback.format_exc()
            self.msgBox(message)
    
    def setPlotButton(self):
        #enables button and combobox to change species plot
        self.changePlotButton.setEnabled(True)
    
        #sets species names in combobo which change plot
        self.spPlotBox.clear()
        self.spPlotBox2.clear()
        self.SpStr2 = self.SpStr.copy()
        self.SpStr2.append('All')
        self.spPlotBox.addItems(self.SpStr2)
        self.spPlotBox2.addItems(self.SpStr)
    
    def computeInput(self):
        #computes input
        uValve, expTime = s2b.pulse_expr(np.array(self.valTon),
                                             np.array(self.valTdur), 
                                             np.array([self.expTdur], float))
        if self.inputT == 'Pulse':
            #takes valve values as input
            inpShock = uValve
            print(inpShock.shape)
        elif self.inputT == 'Model':
            self.inputModelW.regressorInp["inpU"] = uValve
            self.inputModelW.regressorInp["Vtime"] = expTime
            print(self.inputModelW.regressorInp)
            print(uValve.shape, expTime.shape)
            
            #computes system input
            ExprInp = s2b.solveODE(self.inputModelW.regressorInp["ParValues"], self.inputModelW.regressorInp)
            inpShock = ExprInp[-1,:]
        return expTime, inpShock
           
    def plotting(self, expTime, inpShock, AllcellExpr, idxPlot): 
        #sets if there is a secondary plot
        if self.radioPlot2.isChecked():
            self.graphPlot.resize(891,271)
            self.graphPlot3.resize(891,271)
        else:
            self.graphPlot.resize(891,411)
            self.graphPlot3.resize(891,411)
            
        #shows plotting state
        self.labeliter.setText("Plotting...")
        QApplication.processEvents()
        
        #plots input system
        self.graphPlot2.canvas.axes.plot(expTime, inpShock, 'r')
        self.graphPlot2.canvas.axes.fill_between(expTime, inpShock, color='r')
        
        #gets species to be plot
        print(self.outPlot.currentText())
        print(idxPlot)
        
        #plots when there is variability
        if self.noiseBox.isChecked() or self.intrinsicBox.isChecked() or self.extrinsicBox.isChecked() or len(AllcellExpr.shape) == 3:
            if self.indPlotEnable.isChecked():
                for i in range(0, self.nCells):
                    self.graphPlot.canvas.axes.plot(expTime, AllcellExpr[i,:,idxPlot], linewidth=1)
                    
                    #secondary species plot
                    if self.radioPlot2.isChecked() and self.subplotSP == True:
                        self.graphPlot4.canvas.axes.plot(expTime, AllcellExpr[i,:,self.indexPlot2], 
                                                         linewidth=1)  
                    
                    #shows simulation status
                    self.labeliter.setText("Plotting... " + str(i+1) + '/' + str(self.nCells))
                    QApplication.processEvents()
                #end for i
            #end if self.indPlotEnable
            
            #adjusts data to be contain between the 95 % of the output. (2.5-97.5)
            ##########################################################
            ExprMin = np.quantile(AllcellExpr[:,:,idxPlot], 0.025, axis=0)
            ExprMax = np.quantile(AllcellExpr[:,:,idxPlot], 0.975, axis=0)
            ##########################################################
            
            #plots population output
            self.graphPlot3.canvas.axes.plot(expTime, np.median(AllcellExpr[:,:,idxPlot], axis=0), label=self.SpStr[idxPlot], linewidth=2)
            self.graphPlot3.canvas.axes.legend(loc='best')
            self.graphPlot3.canvas.axes.fill_between(expTime, ExprMin,ExprMax, alpha = 0.4)
            
            #highlights negative values in the plot
            if np.amin(AllcellExpr[:,:,idxPlot]) < 0:
                self.graphPlot3.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,idxPlot]),
                                                 color = 'C3', alpha = 0.1)
                self.graphPlot.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,idxPlot]),
                                                 color = 'C3', alpha = 0.1)
            #plots secondary population species
            if self.radioPlot2.isChecked() and self.subplotSP == True:
                #adjusts data to be contain between the 95 % of the output. (2.5-97.5)
                ##########################################################
                ExprMin2 = np.quantile(AllcellExpr[:,:,self.indexPlot2], 0.025, axis=0)
                ExprMax2 = np.quantile(AllcellExpr[:,:,self.indexPlot2], 0.975, axis=0)
                ##########################################################

                #plots population output 
                self.graphPlot5.canvas.axes.plot(expTime, np.median(AllcellExpr[:,:,self.indexPlot2], axis=0), label=self.SpStr[self.indexPlot2], linewidth=2)
                self.graphPlot5.canvas.axes.legend(loc='best')
                self.graphPlot5.canvas.axes.fill_between(expTime, ExprMin2,ExprMax2, alpha = 0.4)                 
                
                #highligths negative values in the plots
                if np.amin(AllcellExpr[:,:,:]) < 0:
                    self.graphPlot4.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,self.indexPlot2]),
                                                     color = 'C3', alpha = 0.1)
                    self.graphPlot5.canvas.axes.fill_between(expTime, np.amin(AllcellExpr[:,:,self.indexPlot2]),
                                                     color = 'C3', alpha = 0.1)
        else:
            #plots when there is not variability
            self.graphPlot.canvas.axes.plot(expTime, AllcellExpr[idxPlot,:], 
                                            label=self.SpStr[idxPlot], linewidth=1)
            self.graphPlot.canvas.axes.legend(loc='best')
            
            #plots population output
            self.graphPlot3.canvas.axes.plot(expTime, AllcellExpr[idxPlot,:], linewidth=2, label=self.SpStr[idxPlot])
            self.graphPlot3.canvas.axes.legend(loc='best')
            
            #secondary species plot
            if self.radioPlot2.isChecked() and self.subplotSP == True:
                self.graphPlot4.canvas.axes.plot(expTime, AllcellExpr[self.indexPlot2,:], 
                                                 label=self.SpStr[self.indexPlot2], linewidth=1)
                self.graphPlot4.canvas.axes.legend(loc='best')
                
                #plots population output
                self.graphPlot5.canvas.axes.plot(expTime, AllcellExpr[self.indexPlot2,:], linewidth=2, label=self.SpStr[self.indexPlot2])
                self.graphPlot5.canvas.axes.legend(loc='best')

        #draws plots into objects
        self.graphPlot.canvas.draw()
        self.graphPlot2.canvas.draw()
        self.graphPlot3.canvas.draw()
        self.graphPlot4.canvas.draw()        
        self.graphPlot5.canvas.draw()        
        
        #disables secondary species plot
        self.subplotSP = False
        
        #shows plotting status
        self.labeliter.setText("Plot Done!")
        QApplication.processEvents()
  
    def computeSim(self, nCells, regressor): 
        #tag to show simulation status
        self.labeliter.setText("Simulating...")
        QApplication.processEvents()

        #array to save individual outputs
        AllcellExpr = np.zeros((nCells, len(regressor["inpU"]), 
                                len(regressor["ODEs"])))

        #system simulation depending on the variability type
        if self.noiseBox.isChecked() or self.intrinsicBox.isChecked() or self.extrinsicBox.isChecked():
            #process extrinsic variability
            if self.extrinsicBox.isChecked():
                print("MxEf")
                #computes natural log from parameter values
                self.logVpars = np.log(regressor["vPars"])

                #creates log-covariance matrix
                self.covStr = self.covMatrix.toPlainText()
                self.covStr = self.covStr.split(';')
                self.covStr = [x.split(',') for x in self.covStr]
                
                #converts logcovariance matrix
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
                self.multipars = np.exp(np.random.multivariate_normal(self.logVpars,
                                                                      self.valCov, 
                                                                      nCells))
                print(self.logVpars)
                print(self.valCov)
                print(self.multipars)
            else:
                #array with original parameters
                self.multipars = np.tile(regressor["vPars"], (nCells,1))
            print(self.multipars)

            t2 = time.time()            
            #population simulation
            for k in range(0,nCells):
                print('Célula: ' + str(k+1))
                #tag to show simulation iteration
                self.labeliter.setText("Simulating... " + str(k+1) + '/' + str(nCells))
                QApplication.processEvents()
            
                #processes intrinsic variability
                if self.intrinsicBox.isChecked():
                    print("Gill")
                    #numeric stochastic solution of differential equations
                    try:
                        cellExpr = s2b.cGill(self.multipars[k,:], regressor)
                    except:
                        cellExpr = s2b.gillAl(self.multipars[k,:], regressor)
                else:
                    #numeric deterministic solution of differential equations
                    try:    
                        t0 = time.time()
                        cellExpr = s2b.solveODEpy(self.multipars[k,:], regressor)
                        t1 = time.time()
                        print("Runs in", t1-t0)
                        
                    except:
                        cellExpr = s2b.solveODE(self.multipars[k,:], regressor)
                    
                #stores individual simulated outputs
                for j in range(0, len(cellExpr)):
                    AllcellExpr[k,:,j] = cellExpr[j,:]
                #end for j
            #end for k
            
            t3 = time.time()
            print("Whole population runs in", t3-t2)
            
            if self.noiseBox.isChecked():
                print("Noise")
                
                #gets noise parameters
                self.changeErrorPar()
                
                #normal distribution
                mu, sigma = 0, 1 #mean and standard deviation
                ndis = np.random.normal(mu, sigma, (nCells,len(regressor["inpU"])))
                
                #computes individual output with noise for each species
                hvar = self.a_err + self.b_err*AllcellExpr[:,:,-1]     
                noisySp = AllcellExpr[:,:,-1] + hvar*ndis
                
                #store species profiles with noise
                AllcellExpr[:,:,-1] = noisySp   
        else:
            #deterministic output
            #numeric solution of ODEs
            try:
                cellExpr = s2b.solveODEpy(regressor["vPars"], regressor)
            except:
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
                #compares matrix values with its equivalent
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
            
        #shows number of cells
        self.outODEs.append("Number of Cells: " + str(nCells))

    #function to update experiment duration value
    def changeTdur(self):
        self.expTdur = self.spinTdur.value()
        
    #gets noise parameters
    def changeErrorPar(self):
        self.a_err = self.spinAval.value()
        self.b_err = self.spinBval.value()
        print(self.a_err, self.b_err)
        
    #function to update number of cells
    def changeCells(self):
        self.nCells = self.numCells.value()
        if self.nCells > 1:
            self.sampleCells.setMaximum(self.nCells)
    
    #opens window to define input model
    def setInputModel(self):
        self.inputModelW.exec_()
    
    #enables input box to define log-covarince matrix
    def enable_cov(self):
        if self.extrinsicBox.isChecked():
            self.covMatrix.setEnabled(True)
        else:
            self.covMatrix.setEnabled(False)
    
    #enables input to define noise parameter values
    def enable_noise(self):
        if self.noiseBox.isChecked():
            self.spinAval.setEnabled(True)
            self.spinBval.setEnabled(True)
        else:
            self.spinAval.setEnabled(False)
            self.spinBval.setEnabled(False)
               
    #changes species to be plotted
    def changePlot(self):
        #clean plots
        self.clearPlot()
        #disable while plitting
        self.changePlotButton.setEnabled(False)
                
        #gets species to be plotted
        print(self.spPlotBox.currentText())
        self.idxPlot2 = self.SpStr2.index(self.spPlotBox.currentText())
        print(self.idxPlot2)
        
        #gets secondary species to be plotted
        if self.radioPlot2.isChecked():
            self.indexPlot2 = self.SpStr.index(self.spPlotBox2.currentText())
            print(self.indexPlot2)
            print(self.spPlotBox2.currentText())
            #enable secondary species subplot
            self.subplotSP = True
        
        #plots chose species
        if not(self.spPlotBox.currentText() == "All"):
            self.plotting(self.expTime, self.inpShock, self.AllcellExpr, self.idxPlot2)
        
        else:
            #plots all outputs
            for i in range(len(self.SpStr)):
                self.plotting(self.expTime, self.inpShock, self.AllcellExpr, i)
        
        #enables after plotting
        self.changePlotButton.setEnabled(True)
        #disables secondary species subplot
        self.subplotSP = False
    
    #enables combobox to choose secondary species to be plotted
    def enablePlot2(self):
        if self.radioPlot2.isChecked():
            self.spPlotBox2.setEnabled(True)
        else:
            self.spPlotBox2.setEnabled(False)
                    
    def inferate(self):
        #restores size of main plot
        self.graphPlot.resize(891,411)
        self.clearPlot()

        #disables while infering
        self.enableSimOpt(False, False)

        #sets tab 1
        self.tabWidget.setCurrentIndex(0)
        
        #change tab names
        self.tabWidget.setTabText(0, "Profile Fitting")
        self.tabWidget.setTabText(1, "Function Value")

        #gets type modelo to perform inference
        self.modelInfr = self.modelBox.currentText()
        print(self.modelInfr)
        
        #takes initial parameter values
        self.iparStr = self.pars0.text()
        self.iparStr = re.findall('[0-9.0-9]+', self.iparStr)
        
        try:
            self.beta0 = list(map(float, self.iparStr))
            
            #check parameter values
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
            
            #computes inferring process
            self.betacal = self.computeInfr(self.modelInfr, self.beta0, self.Obs, self.regressor)
            
            #shows information estimated during inferring process
            self.showInfoInfr(self.betacal)
            
            #save parameter evolution through iteration for each cell
            try:
                if self.parAny.isChecked():
                    self.clearSubplot()
                    self.parAnalysis(self.iterPars, self.nparStr + ['a', 'b'],
                                         np.concatenate((self.valPar, np.array([self.a_err, self.b_err])), axis=None))
            except:
                self.outODEs.append("\nThe number of iterations required to estimate each cell was different.")

            self.labeliter.setText("Inference Done!")
            
            #enable when inference is done
            self.enableSimOpt(True, True)
            
        except:
            self.enableSimOpt(True, True)
            self.outODEs.append("\nPlease, review the input information.")
            message = traceback.format_exc()
            self.msgBox(message)

    #perform model stimations
    def computeInfr(self, model, beta0, Obs, regressor):
        #lets know inferring state
        self.stateN = 0
        self.labeliter.setText("Inferring...")
        QApplication.processEvents()                   
        
        #determine mean and standar deviation curves
        if not np.isnan(Obs).any():
            self.uObs = np.mean(Obs, axis=0)
            self.sdObs = np.std(Obs, axis=0)
        else:
            self.uObs = np.nanmean(Obs, axis=0)
            self.sdObs = np.nanstd(Obs, axis=0)
            
        #initial noise parameter values
        self.A0 = self.A0Box.value()
        self.B0 = self.B0Box.value()
        self.errAB = np.array([self.A0, self.B0])
        print(self.errAB)
        
        #number of iterations to be perform
        self.iterInfr = self.iterBox.value()
        
        #value tolerated between iteration for convergence
        self.tolerance_value = self.spinTolerance.value()
        
        #clear costarray and parsArray
        self.costArray = []
        self.parsArray = []
        
        #get algorithm used to estimate parameters
        self.algorithm = self.AlgoBox.currentText()
        
        if self.algorithm == "CMA-ES":
            try:
                self.getSig = QInputDialog.getDouble(self, 'CMA-ES Parameter', 'Enter Sigma value:', 0.3)
                self.sigma = self.getSig[0]
            except:
                self.sigma = 0.3
            
        self.step = 0
        if model == "Mean Curve Fit":
            #use curve_fit functon to perform regression over the mean curve
            #to infer parameters
            def ODEfit(regressor, *pars):
                #solve simbolic ODEs using new parameters
                cellExpr = s2b.solveODE(pars, regressor)                
                print(pars)
                #iterations count
                self.step += 1
                if self.step % 20:
                    #shows inferring status
                    self.inferState()
                    self.plotInfr(self.expTime, self.uObs, cellExpr[-1,:])
                return cellExpr[-1,:]
            
            betacal, cov = curve_fit(ODEfit, regressor, self.uObs, beta0)
            print(betacal, len(betacal))
        
        elif model == "Mean Cell":
            #uses likelihood function to estimate population parameters
            regressor["Obs"] = Obs
            GuessPars = np.concatenate((beta0, self.errAB))
            
            #likelihood model function
            self.step = 0
            def MLLWrap(pars, regressor):
                nPars = pars[-2:]
                
                expr = s2b.solveODEpy(pars[:-2], regressor)
                uS = expr[-1,:]
                
                sdS = nPars[0] + nPars[1]*uS
                
                fSim = np.tile(uS, (len(regressor["Obs"]),1))
                hSim = np.tile(sdS,(len(regressor["Obs"]),1))
                
                mll = s2b.MLLmeasure(regressor["Obs"], fSim, hSim)
                print(mll)

                #iterations count
                self.step += 1
                if not self.step % 20:
                    #show inferring state
                    self.inferState()
                    self.plotInfr(self.expTime, self.uObs, uS, self.sdObs, sdS)
                    if self.liveInfBox.isChecked():
                        self.costPlot(mll)
                        self.parsPlot(pars)
                return mll
            
            if self.algorithm == "Downhill Simplex":
                minimum, allvecs = fmin(MLLWrap, GuessPars, args=(regressor,), 
                               xtol=self.tolerance_value, 
                               ftol=self.tolerance_value,
                               maxiter=self.iterInfr, retall=True)
                
                niter = len(allvecs)
                
#                res = minimize(MLLWrap, GuessPars, method='Nelder-Mead', 
#                               tol=self.tolerance_value, args=(regressor,), 
#                               options={'maxiter':self.iterInfr,
#                                        'return_all': True})
#                minimum = res.x
#                allvecs = res.allvecs
                
            elif self.algorithm == "CMA-ES":
                
                res = cma.fmin(MLLWrap, GuessPars, self.sigma, args=(regressor,), 
                               options={'maxiter':self.iterInfr,
                                        'tolx': self.tolerance_value, 
                                        'tolfun': self.tolerance_value})
                
#                res = cma.fmin(MLLWrap, GuessPars, self.sigma, args=(regressor,), 
#                               options={'maxiter':self.iterInfr,
#                                        'tolx': self.tolerance_value, 
#                                        'tolfun': self.tolerance_value,
#                                        'bounds': [0, None]})
                minimum = res[0]
                niter = res[4]
            
            if self.liveInfBox.isChecked():
                self.costPlot(0, True, niter)
                self.parsPlot(0, True, niter)
            
            #compute and plot system with parameters found.
            cellCal = s2b.solveODEpy(minimum[:-2], regressor)
            noiseAB = minimum[-2:]
            sdcellCal = noiseAB[0] + noiseAB[1]*cellCal[-1,:]
            self.plotInfr(self.expTime, self.uObs, cellCal[-1,:], self.sdObs, sdcellCal)
            print(minimum)
            
            #get parameters values
            betacal = minimum[:-2]
            self.noiseCal = minimum[-2:]
            self.noiseCal = np.round(self.noiseCal, 5)
            
            #save parameter used in each iteration during estimation
            if self.parAny.isChecked():
                self.iterPars = np.array(allvecs)

        elif model == "Moment Based":
            #use KLD to fit moments and estimate parameters
            GuessPars = np.concatenate((beta0, self.errAB))

            #compute moments
            #determine system without input
            if not self.radioEntrada.isChecked():
                self.ODE2M, self.vars2M = s2b.simbMoments(self.SpStr,
                    self.valReac, self.valProd,self.nparStr)
                
                s2b.model2MDefiner(self.vars2M["nameVar"], self.ODE2M, 
                               self.vars2M["pars"])
            #determine system with input
            else:
                self.ODE2M, self.vars2M = s2b.simbMoments(self.SpStr,
                    self.valReac, self.valProd, self.nparStr, 
                    inputN=self.nameInp, indexU = self.idxReaction)
                
                s2b.model2MDefiner(self.vars2M["nameVar"][1:], self.ODE2M, 
                               self.vars2M["pars"])
            
            print(self.ODE2M)
            
            #moment initial conditions
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
                    "uObs":self.uObs,
                    "sdObs":self.sdObs
                    }
            regressor2.update(self.vars2M)
            
            #find second moment index of target species
            regressor2["idx2M"] = -1
            regressor2["errPars"] = self.errAB
            
            #KLD cost function  
            self.step = 0
            def KLDmomentsWrap(Allpars, regressor):
                uM, sdM = s2b.solve2M(Allpars[:-2], Allpars[-2:], regressor)
                
                mcKLD = s2b.KLDmeasure(regressor["uObs"], regressor["sdObs"], uM, sdM)
                print(mcKLD)
                
                self.step += 1
                if not self.step % 20:
                    #shows inferring status
                    self.inferState()
                    self.plotInfr(self.expTime, self.uObs, uM, self.sdObs, sdM)
                    if self.liveInfBox.isChecked():
                        self.costPlot(mcKLD)
                        self.parsPlot(Allpars)
                return mcKLD
            
            if self.algorithm == "Downhill Simplex":
                minimum, allvecs = fmin(KLDmomentsWrap, GuessPars, 
                                args=(regressor2,),
                                xtol=self.tolerance_value, 
                                ftol=self.tolerance_value,
                                maxiter=self.iterInfr, retall=True)
                niter = len(allvecs)
                
#                res = minimize(KLDmomentsWrap, GuessPars, method='Nelder-Mead', 
#                               tol=self.tolerance_value, args=(regressor2,), 
#                               options={'maxiter':self.iterInfr,
#                                        'return_all': True})
#                minimum = res.x
#                allvecs = res.allvecs
            
            elif self.algorithm == "CMA-ES":
                
                res = cma.fmin(KLDmomentsWrap, GuessPars, self.sigma, 
                               args=(regressor2,), 
                               options={'maxiter':self.iterInfr,
                                        'tolx': self.tolerance_value, 
                                        'tolfun': self.tolerance_value})
                
#                res = cma.fmin(KLDmomentsWrap, GuessPars, self.sigma, 
#                               args=(regressor2,), 
#                               options={'maxiter':self.iterInfr,
#                                        'tolx': self.tolerance_value, 
#                                        'tolfun': self.tolerance_value,
#                                        'bounds': [0, None]})
                minimum = res[0]
                niter = res[4]
            
            if self.liveInfBox.isChecked():
                self.costPlot(0, True, niter)
                self.parsPlot(0, True, niter)
            
            #get parameters inferenced
            betacal  = minimum[:-2]
            print(minimum)
            self.noiseCal = minimum[-2:]
            self.noiseCal = np.round(self.noiseCal, 5)
            
            #plots curve with infered parameters
            uC, sdC = s2b.solve2M(betacal, self.noiseCal, regressor2)
            self.plotInfr(self.expTime, self.uObs, uC, self.sdObs, sdC)
            
            #save parameter used in each iteration during estimation
            if self.parAny.isChecked():
                self.iterPars = np.array(allvecs)
        
        elif model == "Two-Stage":
            #uses likelihood function to estimate population parameters
            GuessPars = np.concatenate((beta0, self.errAB))
            self.samplesC = self.sampleCells.value()
            
            #check sample amount to be bigger than one
            if self.samplesC < 2 and len(Obs) < 2:
                self.outODEs.append("Population is too small.")
                QApplication.processEvents()                   
                return
            
            #subpopulation
#            rIdx = np.random.choice(Obs.shape[0], self.samplesC, replace= False)
            rIdx = np.linspace(0, Obs.shape[0]-1, self.samplesC, dtype=int)
            print(rIdx)
            self.subPop = Obs[rIdx,:]
            print(self.subPop.shape)
            
            #array to store single estimated parameters
            popPars = np.zeros((self.samplesC, len(beta0)))
            popA = np.zeros(self.samplesC)
            popB = np.zeros(self.samplesC)
            
            #likelihood function
            self.step = 0
            def MLLWrap(pars, regressor):
                nPars = pars[-2:]
                
                expr = s2b.solveODEpy(pars[:-2], regressor)
                uS = expr[-1,:]
                
                sdS = nPars[0] + nPars[1]*uS
                
                fSim = np.tile(uS, (len(regressor["Obs"]),1))
                hSim = np.tile(sdS,(len(regressor["Obs"]),1))
                
                mll = s2b.MLLmeasure(regressor["Obs"], fSim, hSim)

                #iterations count
                self.step += 1
                if not self.step % 20:
                    print(mll)
                    #show inferring state
                    self.inferState()
                    self.plotInfr(self.expTime, 
                                  self.subPop[self.iCell,:], 
                                  uS, np.zeros(len(self.expTime)), sdS)
                    if self.liveInfBox.isChecked():
                        self.costPlot(mll)
                        self.parsPlot(pars)
                return mll
            
            #create an array to store parameters if it is required
            if self.parAny.isChecked():
                self.iterPars = np.zeros((self.iterInfr, len(GuessPars), self.samplesC))
            
            #individual estimation of cell population by two-stage
            for self.iCell in range(self.samplesC):
                #reset array of cost values
                self.costArray = []
                
                #change sample to perform estimation
                regressor["Obs"] = np.array([self.subPop[self.iCell,:]])
                
                #estimation function
                if self.algorithm == "Downhill Simplex":
                    minimum, allvecs = fmin(MLLWrap, GuessPars, args=(regressor,),
                                xtol=self.tolerance_value, ftol=self.tolerance_value,
                                maxiter=self.iterInfr, retall=True)
                    niter = len(allvecs)
                    
#                    res = minimize(MLLWrap, GuessPars, method='Nelder-Mead', 
#                               tol=self.tolerance_value, args=(regressor,), 
#                               options={'maxiter':self.iterInfr,
#                                        'return_all': True})
#                    minimum = res.x
#                    allvecs = res.allvecs
                    print(minimum)
                    
                elif self.algorithm == "CMA-ES":
                    
                    res = cma.fmin(MLLWrap, GuessPars, self.sigma, 
                                   args=(regressor,), 
                                   options={'maxiter':self.iterInfr,
                                            'tolx': self.tolerance_value, 
                                            'tolfun': self.tolerance_value})
                    
#                    res = cma.fmin(MLLWrap, GuessPars, self.sigma, 
#                                   args=(regressor,), 
#                                   options={'maxiter':self.iterInfr,
#                                            'tolx': self.tolerance_value, 
#                                            'tolfun': self.tolerance_value,
#                                            'bounds': [0, None]})
                    minimum = res[0]
                    niter = res[4]
            
                if self.liveInfBox.isChecked():
                    self.costPlot(0, True, niter)
                    self.parsPlot(0, True, niter)
                
                #save parameters use in each iteration during estimation of
                #each cell
                try:
                    if self.parAny.isChecked():
                        self.iterPars[:,:,self.iCell] = np.array(allvecs)
                except:
                    #sometimes the number of iterations required to estimate
                    #each cell if diffenrent. Therefore, it will be no possible
                    #to perform parameter analysis
                    pass
                    
                #save estimated parameters
                popPars[self.iCell,:] = minimum[:-2]
                popA[self.iCell] = minimum[-2]
                popB[self.iCell] = minimum[-1]
            #end for iCell
            
            popAB = np.array([popA, popB])
            
            betacal, self.covPars, self.noiseCal = self.tsStats(np.absolute(popPars), 
                                                                np.absolute(popAB))
            self.covPars = np.round(self.covPars, 5)
            print(betacal)
            print(self.covPars)
            print(self.noiseCal)
            
            newMultiPars = np.exp(np.random.multivariate_normal(np.log(betacal), 
                                                                self.covPars, 
                                                                len(Obs)))
            
            ##############################
            popExpr = np.zeros((len(Obs), len(self.expTime)))
            ##############################
            
            for newPars in range(len(Obs)):
                cellCal = s2b.solveODEpy(newMultiPars[newPars,:], regressor)
                popExpr[newPars, :] = cellCal[-1,:]
                
                self.labeliter.setText("Simulating... " + str(newPars + 1) + '/' + str(len(Obs)))
                QApplication.processEvents()
            
            ##############################################################3
            #normal distribution
            mu, sigma = 0, 1 #mean and standard deviation
            ndis = np.random.normal(mu, sigma, (len(Obs),len(self.expTime)))
            
            #computes individual output with noise for each species
            hvar = self.noiseCal[0] + self.noiseCal[1]*popExpr     
            popExpr = popExpr + hvar*ndis
            ############################################################
            
            ucellCal = np.mean(popExpr, axis=0)
            sdcellCal = np.std(popExpr, axis=0)
            
            #a partir de los parametros estimados, simula un individuo, toma
            #la media y aplica modelo de ruido para calcular varianza. Usa 
            #ruido gausiano junto al modelo de ruido.
            self.plotInfr(self.expTime, self.uObs, ucellCal, self.sdObs, 
                          sdcellCal, True, self.subPop)
        
            self.noiseCal = np.round(self.noiseCal, 5)
        return np.round(betacal, 5)
    
    def parAnalysis(self, iterPars, nparStr, ParsRef):
        
        #get the current date and time
        now = datetime.now()
        dt_string = now.strftime("%d_%m_%Y-%H_%M_%S")
        
        #create a folder to save plots
        try:
            os.makedirs('parameter_analysis/')
        except:
            pass
        
        #create a folder to save current analysis
        os.makedirs('parameter_analysis/' + dt_string + '/')
        
        #save parameter in a .npy binary file
        np.save('parameter_analysis/' + dt_string + '/iteration_parameters.npy', iterPars)
        
        #get dimenstions from parameter array
        dimPars = iterPars.shape
        x_axis = np.linspace(1, dimPars[0], dimPars[0])
        
        #create and save plot for each cell
        if len(dimPars) == 3:
            fig, axs = plt.subplots(dimPars[1], sharex=True)
            for i in range(dimPars[2]):
                #create figure for each cell
                fig.suptitle('Cell ' + str(i+1))
                #plot all parameter evolution
                for j in range(dimPars[1]):
                    axs[j].cla()
                    axs[j].plot(x_axis, iterPars[:, j, i])
                    axs[j].grid()
                    axs[j].set(ylabel=nparStr[j])
                plt.xlabel('Iterations')
                
                #save plots as png image
                self.labeliter.setText("Saving... " + str(i + 1) + '/' + str(dimPars[2]))
                QApplication.processEvents()
                plt.savefig('parameter_analysis/' + dt_string + '/cell_' + str(i+1) + '.png')
            #end for i
        
        #plot population parameters iterations
        self.tempAxes = []
        subIdx = str(dimPars[1]) + '1'
        
        if len(dimPars) == 2:
            for k in range(dimPars[1]):
                tempPar = iterPars[:, k]
                tempAx = self.graphPlot6.canvas.figure.add_subplot(int(subIdx + str(k+1)))
                self.tempAxes.append(tempAx)
                tempAx.clear()
                tempAx.plot(x_axis, tempPar, color='tab:red')
                tempAx.plot(np.array([x_axis[0], x_axis[-1]]), 
                   np.array([ParsRef[k], ParsRef[k]]), 'k')
                tempAx.annotate(str(ParsRef[k]), xy=(x_axis[-1],ParsRef[k]), 
                   fontsize=8)
                tempAx.grid()
                tempLab = str(nparStr[k])
                tempAx.set_ylabel(tempLab, fontsize=8)
                tempAx.tick_params(bottom=False, labelbottom=False)
                if k == dimPars[1]-1:
                    tempAx.tick_params(bottom=True, labelbottom=True)
                tempAx.tick_params(axis='both', direction='in', labelsize=8)
            
            self.graphPlot6.canvas.axes.set_xlabel("Iterations", labelpad=16.0, 
                                               fontsize=8)
            self.graphPlot6.canvas.draw()
        
        else:    
            for k in range(dimPars[1]):
                tempPar = iterPars[:, k, :]
                tempAx = self.graphPlot6.canvas.figure.add_subplot(int(subIdx + str(k+1)))
                self.tempAxes.append(tempAx)
                tempAx.clear()
                tempAx.plot(x_axis, np.mean(tempPar, axis=1), color='tab:red')
                tempAx.fill_between(np.linspace(1, dimPars[0], dimPars[0]), 
                   np.mean(tempPar, axis=1) + np.std(tempPar, axis=1),
                   np.mean(tempPar, axis=1) - np.std(tempPar, axis=1), 
                   color='r', alpha = 0.1)
                tempAx.plot(np.array([x_axis[0], x_axis[-1]]), 
                   np.array([ParsRef[k], ParsRef[k]]), 'k')
                tempAx.annotate(str(ParsRef[k]), xy=(x_axis[-1],ParsRef[k]), 
                   fontsize=8)
                tempAx.grid()
                tempLab = str(nparStr[k])
                tempAx.set_ylabel(tempLab, fontsize=8)
                tempAx.tick_params(bottom=False, labelbottom=False)
                if k == dimPars[1]-1:
                    tempAx.tick_params(bottom=True, labelbottom=True)
                tempAx.tick_params(axis='both', direction='in', labelsize=8)
            
            self.graphPlot6.canvas.axes.set_xlabel("Iterations", labelpad=16.0, 
                                               fontsize=8)
            self.graphPlot6.canvas.draw()

        #create and save average parameter evolution
        if len(dimPars) == 2:
            fig, axs = plt.subplots(dimPars[1], sharex=True)
            fig.suptitle("Population Parameter Iteration")
                        
            for k in range(dimPars[1]):
                tempPar = iterPars[:, k]
                axs[k].cla()
                axs[k].plot(x_axis, tempPar, color='tab:red')
                axs[k].plot(np.array([x_axis[0], x_axis[-1]]), 
                   np.array([ParsRef[k], ParsRef[k]]), 'k')
                axs[k].annotate(str(ParsRef[k]), xy=(x_axis[-1],ParsRef[k]), 
                   fontsize=8)
                axs[k].grid()
                axs[k].set(ylabel=nparStr[k])
            #end for k
        else:
            fig.suptitle("Population Parameter Iteration Average and SD")
            
            for k in range(dimPars[1]):
                tempPar = iterPars[:, k, :]
                axs[k].cla()
                axs[k].plot(x_axis, np.mean(tempPar, axis=1), color='tab:red')
                axs[k].fill_between(np.linspace(1, dimPars[0], dimPars[0]), 
                   np.mean(tempPar, axis=1) + np.std(tempPar, axis=1),
                   np.mean(tempPar, axis=1) - np.std(tempPar, axis=1), 
                   color='r', alpha = 0.1)
                axs[k].plot(np.array([x_axis[0], x_axis[-1]]), 
                   np.array([ParsRef[k], ParsRef[k]]), 'k')
                axs[k].annotate(str(ParsRef[k]), xy=(x_axis[-1],ParsRef[k]),
                   fontsize=8)
                axs[k].grid()
                axs[k].set(ylabel=nparStr[k])
        plt.xlabel('Iterations')
        plt.savefig('parameter_analysis/' + dt_string + '/Population_Itearation_Parameter.png')
        
        #Open folder that contains the parameter analysis just created
        nowPath = os.getcwd()
        futPath = nowPath + '/parameter_analysis/' + dt_string
        os.startfile(futPath)
        
    def tsStats(self, parsPop, ABpop):
        #log population parameters
        logPars = np.log(parsPop)

        #mean of log parameters
        mu = np.mean(logPars, axis=0)
        
        #covariance of log parameters
        omega = np.cov(logPars.T)

        #mean of noise parameter        
        ABpars = np.mean(ABpop, axis=1)
        
        return np.exp(mu), omega, ABpars
    
    #lets know that inferring is still ongoing
    def inferState(self):
        if self.stateN == 0:
            self.labeliter.setText("Inferring.")
            
            if self.modelInfr == "Two-Stage":
                self.labeliter.setText("Inferring.   " + str(self.iCell + 1) + "/" + str(self.samplesC))
            
            QApplication.processEvents()
            self.stateN = 1                   

        elif self.stateN == 1:
            self.labeliter.setText("Inferring..")
            
            if self.modelInfr == "Two-Stage":
                self.labeliter.setText("Inferring..  " + str(self.iCell + 1) + "/" + str(self.samplesC))
            
            QApplication.processEvents()                   
            self.stateN = 2
            
        elif self.stateN == 2:
            self.labeliter.setText("Inferring...")
            
            if self.modelInfr == "Two-Stage":
                self.labeliter.setText("Inferring... " + str(self.iCell + 1) + "/" + str(self.samplesC))
            
            QApplication.processEvents()         
            self.stateN = 0
    
    def costPlot(self, costVal, enAxe = False, nAxe = 0):
        if not enAxe == True:
            self.costArray = self.costArray + [costVal]
        else:
            iterAxe = np.linspace(1, nAxe, len(self.costArray))
        
        self.graphPlot3.canvas.axes.clear()
        
        if not enAxe == True:
            self.graphPlot3.canvas.axes.plot(self.costArray)
        else:
            self.graphPlot3.canvas.axes.plot(iterAxe, self.costArray)
            self.graphPlot3.canvas.axes.tick_params(bottom=True, 
                                                    labelbottom=True)
            
        self.graphPlot3.canvas.axes.set_ylabel("Function Value", fontsize=8)
        self.graphPlot3.canvas.axes.grid()        
        self.graphPlot3.canvas.draw()
        QApplication.processEvents()       

    def parsPlot(self, pars, enAxe = False, nAxe = 0):
        self.clearSubplot()
        
        if not enAxe == True:
            #attach new parameters
            self.parsArray = self.parsArray + [pars]   
        else:
            iterAxe = np.linspace(1, nAxe, len(self.parsArray))
        
        #parameters names
        tempName = self.nparStr + ['a', 'b']
        
        #temporal parameters numpy array
        iterPars = np.array(self.parsArray)
        dimPars = iterPars.shape
        
        self.tempAxes = []
        subIdx = str(dimPars[1]) + '1'
#        subIdx = '1' + str(dimPars[1])
        
        for k in range(dimPars[1]):
            tempPar = iterPars[:, k]
            tempAx = self.graphPlot6.canvas.figure.add_subplot(int(subIdx + str(k+1)))
            self.tempAxes.append(tempAx)
            tempAx.clear()
            if not enAxe == True:
                tempAx.plot(tempPar, color='tab:red')
            else:
                tempAx.plot(iterAxe, tempPar, color='tab:red')
                
            tempAx.grid()
            tempLab = str(tempName[k])
            tempAx.set_ylabel(tempLab, fontsize=8)
            tempAx.tick_params(bottom=False, labelbottom=False)
            tempAx.tick_params(axis='both', direction='in', labelsize=8)
        
        if enAxe == True:
            tempAx.tick_params(bottom=True, labelbottom=True)
        self.graphPlot6.canvas.axes.set_xlabel("Iterations", labelpad=16.0, 
                                               fontsize=8)
        self.graphPlot6.canvas.draw()
        QApplication.processEvents()       
    
    #plots infered curve in each iteration
    def plotInfr(self, expTime, uCell, InfrCell, sdCell=np.array([]), sdInfrCell=np.array([]), sampleP = False, sampleSub = np.array([])):
        self.graphPlot.canvas.axes.clear()
        self.graphPlot.canvas.axes.set_ylabel("Concentration", fontsize=8)
        
        #plot observed 
        if self.modelInfr == "Two-Stage":
            self.graphPlot.canvas.axes.plot(expTime, uCell, label="Real")
        else:
            self.graphPlot.canvas.axes.plot(expTime, uCell, label="Observed")
        self.graphPlot.canvas.axes.set_ylim(np.amin(uCell), np.amax(uCell))

        #plotting of output deviation
        if sdCell.size != 0:
            downL = uCell - sdCell
            upL = uCell + sdCell
            self.graphPlot.canvas.axes.plot(expTime, downL, 'b-')
            self.graphPlot.canvas.axes.plot(expTime, upL, 'b-')
            self.graphPlot.canvas.axes.set_ylim(np.amin(downL), np.amax(upL))
        
        #show subpopulation used as sample
        if sampleP == True:
            uSample = np.mean(sampleSub, axis=0)
            sdSample = np.std(sampleSub, axis=0)
            
            self.graphPlot.canvas.axes.plot(expTime, uSample, color='tab:green', label="Observed")
            self.graphPlot.canvas.axes.plot(expTime, uSample + sdSample, 'g')
            self.graphPlot.canvas.axes.plot(expTime, uSample - sdSample, 'g')
        
        #plot inferred population
        self.graphPlot.canvas.axes.plot(expTime, InfrCell, color='orange', label="Inferred")
        
        #plotting of inferred deviation
        if sdCell.size != 0:
            self.graphPlot.canvas.axes.plot(expTime, InfrCell - sdInfrCell, color='tab:orange')
            self.graphPlot.canvas.axes.plot(expTime, InfrCell + sdInfrCell, color='tab:orange')
            
        self.graphPlot.canvas.axes.grid()        
        self.graphPlot.canvas.axes.legend(loc='best')            
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
            
        #show log covariance matrix
        if self.modelInfr == "Two-Stage":
            self.covPars = np.round(self.covPars, 5)
            self.outODEs.append("Log-Covariance Matrix:")
            for i in range(len(self.covPars)):
                    temp = self.covPars[i].tolist()
                    temp = list(map(str, temp))
                    self.outODEs.append('\t|  ' + ', '.join(temp) + '  |')
        QApplication.processEvents()
    
    #function to change inference model and enable specific parameters 
    def changeModel(self):
        if self.modelBox.currentText() == "Mean Curve Fit":
            self.A0Box.setEnabled(False)
            self.B0Box.setEnabled(False)
            self.iterBox.setEnabled(False)
            self.spinTolerance.setEnabled(False)
            self.parAny.setEnabled(False)
            self.AlgoBox.setEnabled(False)
        else:
            self.A0Box.setEnabled(True)
            self.B0Box.setEnabled(True)
            self.iterBox.setEnabled(True)
            self.spinTolerance.setEnabled(True)
            self.parAny.setEnabled(True)
            self.AlgoBox.setEnabled(True)
            
        if  self.modelBox.currentText() == "Two-Stage":
            self.sampleCells.setEnabled(True)
        else:
            self.sampleCells.setEnabled(False)
    
    def changeAlgorithm(self):
        self.algorithm = self.AlgoBox.currentText()
        
        if self.algorithm == "CMA-ES":
            self.parAny.setEnabled(False)
            self.parAny.setChecked(False)
        else: 
            self.parAny.setEnabled(True)
    
    #show error message box
    def msgBox(self, message):
        msg = QMessageBox()
        msg.setWindowTitle("Error!!!")
        msg.setText(message)
        msg.setIcon(QMessageBox.Critical)
        msg.exec_()

    #start and end events of window execution
    def showEvent(self, event):
        self.outODEs.setText("Welcome to SysBioSim!")
    
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