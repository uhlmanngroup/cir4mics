#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:25:03 2023

@author: maria
"""
import yaml
import DeformNPC
import numpy as np
import math

def getVars(config):
    """reads a .yaml file and turns it into a python dictionary
    
    :param config: yaml file to be read 
    :returns: var, a python dictionary version of the config file, which can 
    be modified 
    """
    with open(config) as file:
        var = dict(yaml.load(file, Loader=yaml.FullLoader))
    return var


def getNPCs(var):
    """
    :param var: dictionary of simulation parameters
    returns: Dictionary containing simulated NPCs and their metadata
    """

    return DeformNPC.MultipleNPC(**var)
    # return DeformNPC.MultipleNPC(var["nup"], var["term"], var["model"], n = var["n"], rel = var["rel"], rvar=var["rvar"],
    #                                  thetavar = var["thetavar"], dvar = var["dvar"], symmet = var["symmet"],
    #                                  elliptvar = var["elliptvar"], mag = var["mag"], zmag = var["zmag"], seed = var["seed"])




def getOffsetNPCs(NPCcoords, offset=None, offsetmult=None, justoffset:bool = False):
    """Arrange NPCs on a grid by offsetting them in x and y direction
    :param NPCcoords: Coordinates of NPC, non-offset
    :param offset: distance by which each NPC should be offset. Automatically determined if not provided
    :param offsetmult: multiplier of offset. Default 1.25
    :param justoffset: Only returns offset in nm if True. Default False 
    :returns:
        offsetNPCs: Offset coordinates"""

    offsetmult = 1 if not offsetmult else offsetmult
    offset = offsetmult * np.max(np.abs(NPCcoords[:, :2])) if not offset else offset
  
    if not justoffset:  
        n = int(NPCcoords[-1, 4] + 1)  # number of NPCs
        offsetNPCs = np.copy(NPCcoords)


        # Determine the number of rows and columns needed. The last cells on the grid might stay empty
        ncols = math.ceil(np.sqrt(n))
        nrows = math.ceil(n / ncols)
    
        x = 0  # indexing x coordinate
        y = 1  # indexing y coordinate
        i = 0  # will get updated
    
        for row in range(ncols):
            for col in range(nrows):
                if i < n:
                    offsetNPCs[np.where(offsetNPCs[:, 4] == i), y] += col * 3 * offset
                    offsetNPCs[np.where(offsetNPCs[:, 4] == i), x] += row * 3 * offset
                    i += 1
    
        return offsetNPCs
    return offset


def getNPCcoords(NPCs, var, offset:bool=False, offsetmult = 1, justoffset:bool = False):
    """Generates coordinates of all NPCs with NPCs and nups indexed. 
    :param NPCs: Dictionary containing simulated NPCs and their metadata
    :param var: Dictionary of simulation parameters
    :param offset: NPCs are arranged on a grid if True. Offset can be applied as 
    a separate function if both offset and non-offset NPCs are required. Defaults to False
    :returns: Coordinates of NPCs or coordinates of NPCs and offset coordinates of NPCs of offet = True
    """
    npcs = NPCs["NPCs"]  # Does not yet include tilt and shift
    nupIndex = NPCs["nupIndex"]

    zoffsets = NPCs["zoffsets"]
    NucSideBool = NPCs["nuclear_side_boolean"]

    tiltnucv, tiltcytv = DeformNPC.tiltvectors(var["kappa"], var["n"], var["seed"])
    shiftNuc, shiftCyt = DeformNPC.shiftvectors(
        var["shiftsigma"], var["n"], var["seed"]
    )


    NPCscoords = DeformNPC.MultipleNPCs_coord(
        npcs,
        zoffsets,
        var["symmet"],
        NucSideBool,
        nupIndex,
        tiltnucv,
        tiltcytv,
        shiftNuc,
        shiftCyt,
        seed=var["seed"],
    )
    
    if justoffset:
        return getOffsetNPCs(NPCscoords, justoffset = justoffset)
    
    if offset: 
        return NPCscoords, getOffsetNPCs(NPCscoords, offsetmult=offsetmult)  
    
    return NPCscoords

    
