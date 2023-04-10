#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 10:57:22 2023

@author: maria
"""

from json import dumps
import csv

def printvar(var, indent = 4):
    """ Prints current content of var in an easily readable format 
    
    :param var: dictionary of current model variables 
    :param indent: Indentation level, defaults to 4  
    
    """
    print(dumps(var, indent = indent)) 


def previewCSV(namespace, lines=5):
    count = 0
    with open(namespace) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            print(row)
            count += 1
            if count > lines:
                break