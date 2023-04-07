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

folder = ["shift_mag0/", "shift_mag1/", "shift_mag5/", "shift_mag10/", "shift_mag15/"]
maglist = [0, 1, 5, 10, 15]
export = True 

for i in range(len(folder)):
    data_dir = './dataset/' +folder[i]  # Directory for output files
    
    config = 'config.yaml'
    var = NPC.getVars(config) # Transform config file into a python dictionary 
    
    
    #### Adjust simulation parameters
    var["n"] = 1000# Number of NPCs to be simulated
    NPCi = 0 # index out of n of NPC to be shown in any detail-plot or animation. 0-indexed.  
    var["seed"] = 186571 #seed for reproducibility. Any number but 0
    
    ## Select one or more nups, their terminus, and an NPC model for simulation 
    var["nup"] = ("nup107",)
    var["term"] = ("N",)
    var["model"] = "5a9q" 
    #var["rel"] = True # remove the "#" before var["rel"] = True to select the first nup as reference 
    
    #### Variability parameters
    ### Irregular variability
    var["mag"] = maglist[i] # Magnitude of irregular variability. 0: No deformation, 15: Strong deformation
    var["zmag"] = var["mag"]/2 # magnitude of offset in z. Not computed via springs 
    
    ######################### Geometric variability 
    #var["symmet"] = 8 # Symmetry of simulated NPCs
    
    ## Mean taken from input-model if "None"
    #var["rnew"] = None # Mean radius [nm]. Positive number float or int
    #var["rsigma"] = None # Standard-deviation radius. . Positive number float or int
    
    #var["dnew"] = None # Mean ring distance [nm]. Positive number float or int
    #var["dsigma"] = None # Standard deviation ring distance. Positive number float or int
    
    #var["kappa"] =  None # Controls tilt of individual rings, Kappa of von-mises fisher distribution. Positive number float or int or 0 
    var["shiftsigma"] = 2 # Controls shift of individual rings. Standard deviation of 0-centred normal [nm]
    
    # # Twist between nucleoplasmic and cytoplasmic ring 
    #var["thetanew"] = None # Mean twist angle [rad]
    #var["thetasigma"] = None # standard deviation twist angle
    
    #var["elliptnew"] = None # approximate semiminor/semimajor axis ratio. 
    #var["elliptsigma"] = None # Standard deviation of semiminor/semimajor axis ratio
    
    ########################
    
    #### Run simulations 
    NPCs = NPC.getNPCs(var) # Dictionary of all simulated NPCs
    NPCscoords = NPC.getNPCcoords(NPCs, var) 
    #offset = 1.5 * np.mean(NPCs["rexp"]) # set offset to 1.5 * the expected NPC radius # Not always required 
    #OffsetNPCs = NPC_plotting.OffsetNPCs(NPCscoords, offset) # Not always required 
    
    #### Visualisation 
    # Overview
    #NPC_plotting.plotOverview(OffsetNPCs, NPCs, var, width = 10)
    
    # NPC_plotting.positionVStime(NPCs, index = NPCi, legend = True)
    NPC_plotting.plotDetail(NPCscoords, NPCs, var, index = NPCi, width = 5, mode = "2D")
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
    #NPC_plotting.gethistdata(var, NPCs, featuresAll, featuresElAll, featuresel3DAll, width = 5, bins = None)
    
    #### Export data 
    ## All NPC coordinates and metadata
    if export:
        nameDict, name = exportCSV.names(var)
        exportCSV.MakeCSV(var, NPCs, NPCscoords, nameDict, name, data_dir) 
        
        ## features whole NPC 
        featurescsv = exportCSV.featuresCSV(NPCs, var, name, circle_allrings, ellipse_allrings, data_dir)
        
        ## features per NPC ring 
        exportCSV.featuresCSV_rings(NPCs, var, name, data_dir, circle_allrings, ellipse_allrings)
        
        ## features per subcomplex 
        exportCSV.featuresCSV_subcomplex(NPCs, circle_CRNR, ellipse_CRNR, name, data_dir)