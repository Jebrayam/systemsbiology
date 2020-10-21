# Gene Expression Simulator (GenExpSim) - (ONGOING)
This project was developed in order to help people interested in learn about modeling of gene expression in cell populations. GenExpSim is an informatic tool (desktop application)
created to learn or teach some topics about Systems Biology (Gene Expression). 

## Introduction
Gene expression is a dynamical process through which information encoded in gens is used to  synthesize a genomic product. This product usually is showed as a protein and it is the outcome from a series of biochemical reactions that take place stochastically in a specific biological system. Stochasticity in gen expression affects mainly transcription and translation processes, and it is responsible for the expression levels of a protein differ among a population of genetically identical cells.
The purpous of this project focuses on modeling some of the variability sources that affect how cellular processes response to a stimulus. The variability sources taken into account here come from are extrinsic, intrinsic, and measurement noise. Extrinsic variability arises from phenotypic differences in a cell population. For example, variations in cell size or cell cycle stages. Another factor that affects this type of variability are the cell microenvironment conditions releted to cellular growth, which produce sone adaptative changes in cells. Intrinsic variability takes place as a consequence of the inherent stochastic nature of biochemical reactions associated with transcription, translation, regularization and degradation of mRNA and proteins. Finally, noise measurement take in account how a measuring device is affected by external sources at the time the experiment is performed.

## Software Requirements
The tool is programmed in Python 3.7.4. One of the reasons why it's based on that programming 
language is to make sure that anyone could get access to the tool and if it is the case one could modify or improve it. 

The Python version used to develop this project it is contained into the distribution platform Anaconda, which can be install from its official website. The programming was done mainly using the IDE Spyder 3.3.6, and Qt Designer 5.12.3 was used for the GUI designing. 

To run the tool properly, it is necessary to have installed the following Python packages:
- PyQt5 == 5.12.3
- numpy
- scipy
- sympy
- cython
- matplotlib

## Folders Description
- "GUI" folder contains every file needed to run the tool (application). Make sure each of them are in the same path.

- "Single_Units" folder contains all the algorithms and scripts implemented for each unit. One single unit in this folder can run by itself. The only requirement it is needed is to have the file *simsysbio.py* in the same directory path. In this folder are the functions created and used in this project. Also, a brief explanation about how every script works is given.

## GUI Overview
The *GenExpSim* tool has four main elements: 1) Data input window, 2) System properties window, 3) Command buttons, and 4) System output graphs. The *Data input window* conatains three different sections. Each of these section is meant to execute different functionalities in the tool and each of them depends on the previous section to be executed correctly.

- The **Def** tap refers to the *Define* section. In this section it is possible to define the main properties relate to the biological system. When this section it executed, it shows the differential equations system determined from the data entered by the user and enable the next section.  

- The **Sim** tap refers to the *Simulate* section. There, the user can set and define some simulation settings. The settings include the duration of the experiment, inicial concentration values, number of cell, system input configuration, and define parameters relate to the variability of the output system. When the simulation process finishes, automatically the tool plots the output system. By default it plots firts the last molecular species registered in the previous section, but it is possible to change the species to be plotted and it is necessary one can plot several species at the same time.

- The **Inf** tap refers to the *Infer* Section, which is meant to get the initial values and parameters used to stimate the system output. The variables to be defined are the initial kinetic parameters, the noise parameters (a and b), the number of iterations to be performed, and the model used to run the infering process.

### Main Elements - Define Section 
![def](https://user-images.githubusercontent.com/57733110/96326504-c8ca3a80-0ff6-11eb-8a09-e2568b5728ea.png)

### Simulate Section
![sim](https://user-images.githubusercontent.com/57733110/96326513-dbdd0a80-0ff6-11eb-8cf4-0d4abbe4c2cf.png)

### Infer Section
![inf](https://user-images.githubusercontent.com/57733110/96326516-e13a5500-0ff6-11eb-9682-cbbbfbadf959.png)






