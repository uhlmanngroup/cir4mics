#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 10:57:22 2023

@author: maria
"""

from json import dumps

def printvar(var, indent = 4):
    """ Prints current content of var in an easily readable format 
    
    :param var: dictionary of current model variables 
    :param indent: Indentation level, defaults to 4  
    
    """
    print(dumps(var, indent = indent)) 
