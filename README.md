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

## Files Description
- "Single_Units" folder contains all the algorithms and scripts implemented for each unit. One single unit in this folder can run by itself. The only requirement it is needed is to have the file *simsysbio.py* in the same directory path. In this folder are the functions created and used in this project. Also, a brief explanation about how every script works is given.
- "GUI" folder contains every file needed to run the tool (application). Make sure each of them are in the same path.

## GUI Overview
The *GenExpSim* tool have four main elements: 1) Information input window, 2) System properties window, 3) Command buttons, and 4) System output graphs.
The **Def** tap refers to the *Define* section. In this section it is possible to define the main properties relate to the biological system to be simulated.

![def](https://user-images.githubusercontent.com/57733110/96006984-ab814a80-0e03-11eb-95de-3a3c8ff3311d.png)


