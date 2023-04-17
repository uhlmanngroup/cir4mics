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

Change to the directory folder:  
`cd cir4mics`  

### Option A: Installation from pypi 
`pip install cir4mics` 


### Option B: Installation from Github 
Skip to "Tutorials" if Installation via Option A was successfull  

#### Install poetry 
`pip install poetry`

or following: 
https://python-poetry.org/docs/ 

#### Usage
From the same directory as before 

#### Open a poetry shell
Create or start a virtual environment. Packages installed here will be isolated from the rest of the system, avoiding unfavourable crosstalk.  
`poetry shell`  
The name of the virtual environment should now be indicated in round brackets to the left of the command line prompt. 
If there are issues with this, try first restarting the terminal and then closing any existing environments, such as conda environments.  

#### Install the dependencies 
Dependencies are installed according to pyproject.toml and poetry.lock  
`poetry install`  



#### Run without IDE 
If you are already aware of what's in call_functions.py (See Tutorials)

`poetry run python cir4mics/call_functions.py`  

#### To exit the poetry shell 
`exit` 


## Tutorials 
Run from within the poetry shell and the right working directory (see above).  
### How to simulate NPCs 
`jupyter-lab notebooks/NPC_testlab.ipynb`

#### Run a minimal script  
- Start an IDE of your choice 
- Load cir4mics/scripts/call_functions.py

### How to include new models to the simulator 
`jupyter-lab notebooks/Include_models.ipynb`

### How to simulate temporal dynamics 
Dynamics do not reflect true NPC dynamics, but might be useful when hard to predict dynamics are of interest 

`jupyter-lab notebooks/Dynamics.ipynb` 


## Troubleshoothing (Installation via Option B)
Windows: Git Bash struggles to open poetry shell  
Solution: Use Powershell, cmd.exe, or WSL  
see: https://github.com/python-poetry/poetry/issues/6495  

Windows: `poetry install` is unable to install debugpy  
Solution: run `poetry config installer.modern-installation false`  
see: https://github.com/microsoft/debugpy/issues/1246  

Windows: powershell does not execute `poetry shell` and returns 
`<filepath> cannot be loaded because running scripts is disabled on this system error`  
Solution:  
1. Open the powershell as administrator  
2. Run `Set-ExecutionPolicy RemoteSigned`  
See: https://www.sharepointdiary.com/2014/03/fix-for-powershell-script-cannot-be-loaded-because-running-scripts-is-disabled-on-this-system.html  

Windows: Issues with Microsoft C++ Build Tools (https://visualstudio.microsoft.com/visual-cpp-build-tools/):  
Solution: Use at most version 3.10  

Windows: Powershell upon running `poetry shell` returns `Virtual environment already activated`. 
Solution: Restart the computer 
See: https://visualstudio.microsoft.com/visual-cpp-build-tools/ 


## Contact 
theiss@ebi.ac.uk 
