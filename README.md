## Requirements    
- python, pip, git
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
`git clone https://github.com/MariaTheiss/cir4mics.git`  

Change to the directory folder:  
`cd cir4mics`  

### Install poetry 
`pip install poetry`

or following: 
https://python-poetry.org/docs/ 

## Usage
From the same directory as before 

### Open a poetry shell
Create or start a virtual environment. Packages installed here will be isolated from the rest of the system, avoiding unfavourable crosstalk.  
`poetry shell`  
The name of the virtual environment should now be indicated in round brackets to the left of the command line prompt. 

### Install dependencies 
Dependencies are installed according to pyproject.toml and poetry.lock  
`poetry install`  

### Run the tutorial 
`jupyter-lab cir4mics/NPC_testlab.ipynb`

Everything is properly configured if the tutorial runs to the end.  
The tutorial can be modified for any specific task.  

### Run minimal script 

- Start an IDE of your choice 
- Load cir4mics/cir4mics/call_functions.py

### Run without IDE 
If you are already aware of what's in call_functions.py, and if it doesn't require modification.  

`poetry run python cir4mics/call_functions.py`  


### To exit poetry shell 
`exit` 

## Troubleshoothing 

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


## Tutorials 
### Simulating NPCs: 
`jupyter-lab cir4mics/NPC_testlab.ipynb`

### Including new models into the simulator 
`jupyter-lab cir4mics/Include_models.ipynb`

### Dynamics 
Dynamics do not reflect true NPC dynamics, but might be useful when hard to predict dynamics are of interest 

`jupyter-lab cir4mics/Dynamics.ipynb` 

## Contact 
theiss@ebi.ac.uk 
