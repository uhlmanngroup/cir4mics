#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:25:03 2023

@author: maria
"""
import yaml
import DeformNPC


def getVars(config):
    with open(config) as file: 
        var = dict(yaml.load(file, Loader=yaml.FullLoader))

    return var

def getNPCs(var):
    
    return DeformNPC.MultipleNPC(var["nup"], var["term"], var["model"], n = var["n"], relative = var["rel"], rvar=var["rvar"], 
                                     thetavar = var["thetavar"], dvar = var["dvar"], symmet = var["symmet"], 
                                     elliptvar = var["elliptvar"], mag = var["mag"], zmag = var["zmag"], seed = var["seed"])
    


def getNPCcoords(NPCstest, var):
    
    NPCs = NPCstest["NPCs"] # Does not yet include tilt and shift
    nupIndex = NPCstest["nupIndex"]

    zoffsets = NPCstest["zoffsets"]
    NucSideBool = NPCstest["nuclear_side_boolean"]
    
    tiltnucv, tiltcytv = DeformNPC.tiltvectors(var["kappa"], var["n"], var["seed"])
    shiftNuc, shiftCyt = DeformNPC.shiftvectors(var["shiftsigma"], var["n"], var["seed"]) 
    
    return DeformNPC.MultipleNPCs_coord(NPCs, zoffsets, var["symmet"], NucSideBool, nupIndex, tiltnucv, tiltcytv, shiftNuc, shiftCyt, seed = var["seed"]) 
