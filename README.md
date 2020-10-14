# Gene Expression Simulator (GenExpSim)
This project was developed in order to help people interested in learn about modeling of gene expression in cell populations. GenExpSim is an informatic tool (desktop application)
created to learn or teach some topics about Systems Biology (Gene Expression). 

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

## Folder Description
- "Single_Units" folder contains all the algorithms and scripts implemented for each unit. One single unit in this folder can run by itself. The only requirement it is needed is to have the file *simsysbio.py* in the same directory path. In this folder are the functions created and used in this project. Also, a brief explanation about how every script works is given.
- "GUI" folder contains every file needed to run the tool (application). Make sure each of them are in the same path.

## GUI Overview
The *GenExpSim* tool have four main elements: 1) Data input window, 2) System properties window, 3) Command buttons, and 4) System output graphs. The *Data input window* conatains three different sections. Each of these section is meant to execute different function in the tool and each of them depends on the previous section to be executed correctly.

- The **Def** tap refers to the *Define* section. In this section it is possible to define the main properties relate to the biological system. When this section it executed, it shows the differential equations system determined from the data entered by the user and enable the next section.  
- The **Sim** tap refers to the *Simulate* section. There, the user can set and define some simulation settings. The settings include the duration of the experiment, inicial concentration values, number of cell, system input configuration, and define parameters relate to the variability of the output system. When the simulation process finishes, automatically the tool plots the output system. By default it plots firts the last molecular species registered in the previous section, but it is possible to cheange the species to be plotted and it is necessary one can plot several species at the same time.
- The **Inf** tap refers to the *Infer* Section, which is meant to get the initial values and parameters used to stimate the system output. The variables to be defined are the initial kinetic parameters, the noise parameters (a and b), the number of iterations to be performed, and the model used to run the infering process.

### Main Elements - Define Section 
![def](https://user-images.githubusercontent.com/57733110/96006984-ab814a80-0e03-11eb-95de-3a3c8ff3311d.png)

### Simulate Section
![sim](https://user-images.githubusercontent.com/57733110/96008002-b4bee700-0e04-11eb-9cdd-a4f9126425e4.png)

### Infer Section
![inf](https://user-images.githubusercontent.com/57733110/96008182-e768df80-0e04-11eb-8226-1f159ee523d6.png)






