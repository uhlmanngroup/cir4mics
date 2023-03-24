## Setup 
Requirements: 
- python, pip 
https://pip.pypa.io/en/stable/installation/ 

### Linux, macOS
Open a terminal where you want to install the project  

### Windows 
Recommended Windows command line shells are Windows powershell or cmd.exe. 
Alternatively use Windows Subsystem for Linux (WSL). 


### Intallation 
Clone this repository:  
Tip: Copy/paste commands to avoid typos 
`git clone https://github.com/MariaTheiss/NPCpy.git`  

Change to the directory folder:  
`cd NPCpy`  

### Install poetry 
`pip install poetry`

or following: 
https://python-poetry.org/docs/ 

### Open a poetry shell
Create or start a virtual environment. Packages installed here will be isolated from the rest of the system, avoiding unfavourable crosstalk.  
`poetry shell`
The name of the virtual environment should now be indicated in round brackets to the left of the command line prompt. 

### Install dependencies 
Dependencies are installed according to pyproject.toml and poetry.lock  
`poetry install` 

### Run the tutorial 
`jupyter-lab NPCpy/NPC_testlab.ipynb`

Everything is properly configured if the tutorial runs to the end.  
The tutorial can be modified for any specific task.  

### Run minimal script 

- Start an IDE of your choice 
- Load NPCpy/NPCpy/call_functions.py

### Run without IDE 

`poetry run python NPCpy/call_functions` 


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

Windows: Powershell, upon running `poetry shell` returns `Virtual environment already activated`. 
Solution: Restart the computer 
See: https://visualstudio.microsoft.com/visual-cpp-build-tools/ 


