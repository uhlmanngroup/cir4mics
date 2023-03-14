#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 13:30:07 2023

@author: maria
"""


#working_dir = '/home/maria/Documents/NPCPython' # Working directory (Location of NPC_overview_and_CSV.py)
data_dir = '/home/maria/Documents/NPCPython/deletemedata2/' # Directory for output files
config = 'config.yaml'

#working= os.environ.get("WORKING_DIRECTORY", working_dir)
#if len(sys.argv) > 1: working = sys.argv[1]
#os.chdir( working )

import exportCSV
import NPC 
import NPC_plotting



#var = NPC.variables(working_dir, config)
var = NPC.variables2(config)


var["n"] = 9 # Number of NPCs to be simulated
NPCi = 0 # index of NPC to be shown in any detail-plot or animation 
var["seed"] = 51234 #seed for reproducibility. Any number but 0


var["nup"] = ("nup107",)#, "Nup107", "Nup155")
var["term"] = ("C",)#, "N", "N") #TODO: Force proper syntax. Issues without trailing comma and rel = True
var["model"] = "7r5k"


# ############################
# # Variability parameters
var["mag"] = 15 # Magnitude of irregular variability. 0: No deformation, 15: Strong deformation
var["zmag"] = var["mag"]/2 # magnitude of offset in z. Not computed via springs 


if type(var["nup"]) == str: var["nup"] = tuple([var["nup"]])
if type(var["term"]) == str: var["term"] = tuple(["term"])


 
# functions above not shown in jupyter_lab, but bundling functions from elsewhere in the code. 
# function calls below should go to jupyter lab. 
# So the user should have one function to call for each relevant operation: e.g. generating NPCs, 
# Generating a specific plot, exporting a specific file, etc. 
# Variables can be loaded as-is from config file, or modified. Several config files 
# can be used for different kinds of applications 



NPCstest = NPC.getNPCs(var) # 
NPCscoords = NPC.getNPCcoords(NPCstest, var) 

OffsetNPCs = NPC_plotting.OffsetNPCs(NPCscoords, 80)

NPC_plotting.plotOverview(OffsetNPCs, NPCstest, var)

NPC_plotting.positionVStime(NPCstest, var, 0)
NPC_plotting.plotDetail(NPCscoords, NPCstest, var, NPCi)



NPC_plotting.AnimateDetail(NPCstest, var, NPCi)

NPC_plotting.AnimateOverview(NPCstest, OffsetNPCs, var)


#NPC_overview_and_CSV.MakeCSV(nameDict, name, offset = False) 



nameDict, name = exportCSV.names(var)
exportCSV.MakeCSV(var, NPCstest, NPCscoords, nameDict, name, data_dir) 
