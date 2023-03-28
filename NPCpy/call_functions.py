#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 13:30:07 2023

@author: maria
"""
#### Simulating structurally variable NPCs for microscpy #### 
### 
import exportCSV
import NPC 
import NPC_plotting
import Analyse_deformed
import numpy as np

data_dir = '/home/maria/Documents/NPCPython/deletemedata2/' # Directory for output files

config = 'config.yaml'
var = NPC.getVars(config) # Transform config file into a python dictionary 


#### Adjust simulation parameters
var["n"] = 1 # Number of NPCs to be simulated
NPCi = 0 # index out of n of NPC to be shown in any detail-plot or animation. 0-indexed.  
var["seed"] = 123 #seed for reproducibility. Any number but 0

## Select one or more nups, their terminus, and an NPC model for simulation 
var["nup"] = ("nup107",)
var["term"] = ("C",)
var["model"] = "5a9q"#"simple"
var["model"] = "simple"
#var["rel"] = True # remove the "#" before var["rel"] = True to select the first nup as reference 

#### Variability parameters
### Irregular variability
var["mag"] = 20 # Magnitude of irregular variability. 0: No deformation, 15: Strong deformation
var["zmag"] = var["mag"]/2 # magnitude of offset in z. Not computed via springs 

######################### Geometric variability 
# symmet = 8 # Symmetry of simulated NPCs

## Mean taken from input-model if "None"
# rnew = None # Mean radius [nm]. Positive number float or int
# rsigma = None # Standard-deviation radius. . Positive number float or int

# dnew = None # Mean ring distance [nm]. Positive number float or int
# dsigma = None # Standard deviation ring distance. Positive number float or int

# kappa =  None # Controls tilt of individual rings, Kappa of von-mises fisher distribution. Positive number float or int or 0 
# shiftsigma = None # Controls shift of individual rings. Standard deviation of 0-centred normal [nm]

# # Twist between nucleoplasmic and cytoplasmic ring 
# thetanew = None # Mean twist angle [rad]
# thetasigma = None # standard deviation twist angle

# elliptnew = None # approximate semiminor/semimajor axis ratio. 
# elliptsigma = None # Standard deviation of semiminor/semimajor axis ratio

########################

#### Run simulations 
NPCs = NPC.getNPCs(var) # Dictionary of all simulated NPCs
NPCscoords = NPC.getNPCcoords(NPCs, var) 
offset = 1.5 * np.mean(NPCs["rexp"]) # set offset to 1.5 * the expected NPC radius 
OffsetNPCs = NPC_plotting.OffsetNPCs(NPCscoords, offset)

#### Visualisation 
NPC_plotting.plotOverview(OffsetNPCs, NPCs, var, width = 10)
# NPC_plotting.positionVStime(NPCs, index = NPCi, legend = True)
# NPC_plotting.plotDetail(NPCscoords, NPCs, var, index = NPCi, width = 5, mode = "3D")
# %matplotlib qt # run to view animation. Default is %matplotlib inline
# name = None # 
# NPC_plotting.AnimateOverview(NPCs, OffsetNPCs, var, width = 6, directory = data_dir, name = name, ext = ".gif")

#### Analyse generated NPCs
ellipse_allrings = Analyse_deformed.Ellipses(NPCscoords, membership = NPCs["z_i_all"])
ellipse_CRNR = Analyse_deformed.Ellipses(NPCscoords, membership = NPCs["ringmemall"])

circle_allrings = Analyse_deformed.Circles(NPCscoords, membership = NPCs["z_i_all"])
circle_CRNR = Analyse_deformed.Circles(NPCscoords, membership = NPCs["ringmemall"])  

# Compute further features 
featuresAll, _, _ = Analyse_deformed.meanfeaturesC(NPCs, var, circle_allrings) # Circle features
featuresElAll, _, _ = Analyse_deformed.meanfeaturesE(NPCs, var, ellipse_allrings)
_, _, _, _, featuresel3DAll = exportCSV.col_features(NPCs, circle_CRNR, ellipse_CRNR)

# featuresAll: r, sqsum, residual, zsqsum of circle per NPC 
#featuresel3DAll: centre subcomplex, 1 centre subcomplex 2, centre subcomplex n, tilt subcomplex 1, tilt subcomplex 2 ...

## Show histogram of features 
# NPC_plotting.gethistdata(var, NPCs, featuresAll, featuresElAll, featuresel3DAll, width = 5, bins = None)

#### Export data 
## All NPC coordinates and metadata
nameDict, name = exportCSV.names(var)
exportCSV.MakeCSV(var, NPCs, NPCscoords, nameDict, name, data_dir) 

## features whole NPC 
featurescsv = exportCSV.featuresCSV(nameDict, NPCs, var, name, circle_allrings, ellipse_allrings, data_dir)

## features per NPC ring 
exportCSV.featuresCSV_rings(nameDict, NPCs, var, name, data_dir, circle_allrings, ellipse_allrings)

## features per subcomplex 
exportCSV.featuresCSV_subcomplex(NPCs, nameDict, circle_CRNR, ellipse_CRNR, name, data_dir)