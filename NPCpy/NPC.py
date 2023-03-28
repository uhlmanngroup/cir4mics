#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:25:03 2023

@author: maria
"""
import yaml
import DeformNPC


def getVars(config):
    """reads a .yaml file and turns it into a python dictionary"""
    with open(config) as file: 
        var = dict(yaml.load(file, Loader=yaml.FullLoader))

    return var

def getNPCs(var):
    """Input: var (dictionary)
    Output: Dictionary containing simulated NPCs and their metadata
    """
    return DeformNPC.MultipleNPC(**var)
    # return DeformNPC.MultipleNPC(var["nup"], var["term"], var["model"], n = var["n"], rel = var["rel"], rvar=var["rvar"], 
    #                                  thetavar = var["thetavar"], dvar = var["dvar"], symmet = var["symmet"], 
    #                                  elliptvar = var["elliptvar"], mag = var["mag"], zmag = var["zmag"], seed = var["seed"])
    

def getNPCcoords(NPCs, var):
    """Input: NPCs, var
    """
    npcs = NPCs["NPCs"] # Does not yet include tilt and shift
    nupIndex = NPCs["nupIndex"]

    zoffsets = NPCs["zoffsets"]
    NucSideBool = NPCs["nuclear_side_boolean"]
    
    tiltnucv, tiltcytv = DeformNPC.tiltvectors(var["kappa"], var["n"], var["seed"])
    shiftNuc, shiftCyt = DeformNPC.shiftvectors(var["shiftsigma"], var["n"], var["seed"]) 
    
    return DeformNPC.MultipleNPCs_coord(npcs, zoffsets, var["symmet"], NucSideBool, nupIndex, tiltnucv, tiltcytv, shiftNuc, shiftCyt, seed = var["seed"]) 
