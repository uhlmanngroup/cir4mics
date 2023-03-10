## Setup 

### Required: 

git, conda, a command-line inter interface 

For Windows: 
https://discuss.codecademy.com/t/setting-up-conda-in-git-bash/534473

### Clone repository 

Initialise git in your target folder. 

Next, clone the repository: 

	git clone https://github.com/MariaTheiss/NPCSpringModel.git

### Import environment 

	conda env create -n NPC --file ENV.yml 
	conda activate NPC 

### Load the script 

Run an IDE (e.g. spyder) in the folder that contains NPC_overview_and_CSV.py.  

To open spyder: 

	spyder

Open NPC_overview_and_CSV.py in the IDE

### Adapt NPC_overview_and_CSV.py 

Modify the python code as follows: 

Add the directory that contains NPC-overview_and_CSV.py 

	working dir = "path/to/directory"

Add the path to where output files should go 

	data_dir = "path/to/data/dir" 

## Simulate NPCs

Change the seed to any number except 0 for reproducible results.  

Select the Nup, N or C terminus, and model of interest: 

	nup = "NupOfInterest"
	term = "N" # or "C"
	model = "PDBModelOfInterest"

For a list of available options see preprint. 

Set numbers of NPC to be simulated, 

	n_input = # any integer number

Adjust variability parameters 
Irregular variability: 

	mag = # any number >= 0 

Geometric variability 
Change radius, ring distance, twist angle, and ellipticity here. 
For example: 

	rnew = # change to new mean radius 
	rsigma = # change to standard deviation on a Normal distribution to randomly change the radius 

Note on ellipticity: Due to spring interactions, the input value might not correspond to the output value. 
Generate elliptical NPCs with mag = 0 and read out minor/major axis ratio from an Overview plot with fitted ellipse. 

### Plotting 
Change plotting parameters as follows to not show plots: 

	Overviewplot = {"plot": False , "ellipse": True, "circle": False}
	Detailplot2D = {"plot": False, "showforce" : False}
	XYvsTime = False
	Detailplot3D = False

Overviewplot shows all simulated NPCs with fitted ellipse and/or circle 
Detailplot2D shows the first simulated NPC by default
XYvsTime shows the change in x and y coordinates of nodes of a selected NPC over time. 
This can be used to ensure deformed NPCs are at equilibrium. 
Detailplot3D shows a 3D plot of a selected NPC. 

### Export data and features 

Set to True to export simualted NPCs as a CSV file in the data directory: 

	MakeCSV = False

Set to True to export features of simulated NPCs: 

	featuresCSV = False 

Set to True to export features of individual NPC rings: 

	featuresCSV2 = False
