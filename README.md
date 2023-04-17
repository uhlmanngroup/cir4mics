# CIR4MICS: **c**onfigurable **i**rregular **r**ings **for** **mic**roscopy **s**imulations 
CIR4MICS ("ceramics") is a pipeline that creates artificial datasets of rotationally symmetric structures. 
Out of the box, it generates structurally variable synthetic Nuclear Pore Complexes (NPCs) based on architectural models of the true NPC. 
Users can select one or more N- or C-terminally tagged NPC proteins, and simulate a wide range of geometric variations. 
Rotationally symmetric structures such as the NPC are also represented as a spring-model such that arbitrary deforming forces, 
of user-defined magnitudes, simulate irregularly shaped variations. 

## Requirements    
- python, pip, git, jupyter-lab
https://pip.pypa.io/en/stable/installation/ 

### Linux, macOS
Open a terminal where you want to install the project  

### Windows 
Recommended Windows command line shells are Windows powershell or cmd.exe. 
Alternatively use Windows Subsystem for Linux (WSL). 


## Intallation 
Open a terminal of choice  
Clone this repository:  
Tip: Copy/paste commands to avoid typos  
`git clone https://github.com/uhlmanngroup/cir4mics.git`  

## Install cir4mics 
`pip install cir4mics`

## Run a minimal script 

- Start an IDE of your choice 
- Load cir4mics/cir4mics/call_functions.py

## Tutorials 
Run from within the poetry shell and the right working directory (see above).  
### How to simulate NPCs 
`jupyter-lab cir4mics/NPC_testlab.ipynb`

### How to include new models to the simulator 
`jupyter-lab cir4mics/Include_models.ipynb`

### How to simulate temporal dynamics 
Dynamics do not reflect true NPC dynamics, but might be useful when hard to predict dynamics are of interest 

`jupyter-lab cir4mics/Dynamics.ipynb` 


## Contact 
theiss@ebi.ac.uk 
