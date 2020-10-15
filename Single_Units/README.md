# Single Units Implementation

In this folder are all functions implemented for the project. Here, it is possible to download each function to run it individually. Each function has its own brief explanation and comments related to implementation. 
Some of the functions depend on others to run correctly, so the *simsysbio.py* file is required to be in the same directory path. 

The implemented functions to this project are the next:
- **simbODE**. Determines simbolically the set of diferrential equations that describe the system output.
- **HOGexpr**. This function compute the system input. The output of this function are three signals or profiles. Step signal(Valve), delay step signal (chamber), and model signal (hog). Model signal is based on a osmorregulation process in yeast. This process has its own kinetic parameters, however these can be change to model other processes.
- **solveODE**. Compute the differential equations system. It takes as input the outcome from the *simbODE* function.
