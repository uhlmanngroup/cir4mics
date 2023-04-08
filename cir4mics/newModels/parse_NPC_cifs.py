#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:16:01 2022

@author: maria
"""
import sys
import os
import numpy as np
from gemmi import cif
import pandas as pd
from collections import Counter

# set working directory of input-files and where output files are to be stored
working_dir = "/CoordinatesfromPDB"
working = os.environ.get("WORKING_DIRECTORY", working_dir)
if len(sys.argv) > 1:
    working = sys.argv[1]
os.chdir(working)

# Document of interest
doc = cif.read("7r5k-assembly1.cif")
outputName = "7r5ktest.txt"
output = False
# in the order as they appear in the document (see variable "nups" for reference)
nupnames = [
    "RANBP2",
    "NUP210",
    "ALADIN",
    "NUP93",
    "NUP188",
    "NUP205",
    "NUP155",
    "NDL1",
    "NUP35",
    "NUP54",
    "NUP58NUP45",
    "NUP62",
    "NUP133",
    "NUP107",
    "NUP98NUP96v1",
    "SEC13",
    "SEH1",
    "NUP85",
    "NUP43",
    "NUP160",
    "NUP37",
    "ELYS",
    "NUP98NUP96v2",
    "NUP214",
    "NUP88",
]


block = doc.sole_block()

# Extract relevant data from document
atomID = list(
    block.find_loop("_atom_site.auth_atom_id")
)  # type of atom, we'll later filter for C alpha
seqID = list(
    block.find_loop("_atom_site.auth_seq_id")
)  # index of amino acid of interest
asymID = list(block.find_loop("_atom_site.auth_asym_id"))  # ID of Nup
col_x = list(block.find_values("_atom_site.Cartn_x"))  # Cartesian coordinates of atoms
col_y = list(block.find_values("_atom_site.Cartn_y"))
col_z = list(block.find_values("_atom_site.Cartn_z"))

allvalues = pd.DataFrame([atomID, seqID, asymID, col_x, col_y, col_z]).transpose()
allvalues.columns = ["atomID", "seqID", "asymID", "x", "y", "z"]

allvalues = allvalues[allvalues.atomID == "CA"]  # Filter for only C alpha atoms

# "-" only to be found in Nups IDs of asymmetric unit 2-8, we're interested in asymmetric unit 1
allvalues = allvalues[allvalues.asymID.str.contains("-") == False]

# allvalues = allvalues[allvalues.asymID.str.contains("-5") == True] #opposite rot unit, only to calculate centre

unique_id = np.unique(allvalues.asymID)  # Unique array list of Nup IDs of interest

# Extract indices of N and C terminal resolved amino acids of of C alpha atoms of Nups of interest
indicesC = []
indicesN = []

for value in unique_id:
    indicesN.append(
        allvalues.asymID[allvalues.asymID == value].index[0]
    )  # first aa [0] is most N-terminal
for value in unique_id:
    indicesC.append(
        allvalues.asymID[allvalues.asymID == value].index[-1]
    )  # last aa [-1] is most C-terminal

valuesN = allvalues.filter(items=indicesN, axis=0)
valuesC = allvalues.filter(items=indicesC, axis=0)

nups = np.unique([str(nup[0]) for nup in list(unique_id)])  # nup IDs (types only)
nupsall = np.array(
    [str(nup[0]) for nup in list(unique_id)]
)  # nup IDs, repeated by their frequency

nupsdict = {
    nups[i]: nupnames[i] for i in range(len(nups))
}  # dictionary, keys are nups IDs, entries indicate their names
count = Counter(
    nupsall
)  # dictionary, keys are nup IDs, entries indicate their frequency

if (
    output
):  # print Nup coordinates into a .txt file in a format that can be copy/pasted to Nups_Info.py
    with open(outputName, "w") as f:
        k = 1  # count for nup of a specific type
        ref = ""
        mode = "N"
        idlist = []

        def writeauths(i, values):
            auth = values[values.asymID == i][["asymID"]]
            auth = auth.to_string(index=False, header=False)

            entry = values[values.asymID == i][["x", "y", "z"]]
            entry = entry.to_csv(index=False, header=False).strip()
            return auth, entry

        for i in unique_id:  # for each type of Nup
            if k == 1:  # before a new type of Nup, any terminus (C or N)
                f.write('if (nup == "' + nupsdict[i[0]] + '"):\n')

            #### N-terminus
            if mode == "N":
                idlist.append(i)  # generate list of IDs to later repeat for C-terminus
                auth, entry = writeauths(i, valuesN)
                authnameN = "authN_" + auth

                if k == 1:  # before first Nup of a type, N-terminus
                    f.write('\tif (term == "N"):\n')

                f.write(
                    "\t\t" + authnameN + " = np.array([" + entry + "])\n"
                )  # coordinates

                ref = ref + authnameN + ", "

                # after last Nup of a type, N-terminus
                if k == count[i[0]]:  # counts number of Nups of a specific type
                    f.write(
                        "\t\tref = np.array([" + ref[:-2] + "])\n\n"
                    )  # -2 to remove trailing comma and whitespace
                    k = 0
                    ref = ""
                    mode = "C"  # switch to C terminus once done with N-terminus
                k += 1

            #### C-terminus
            if mode == "C":
                for (
                    i
                ) in (
                    idlist
                ):  # repeat what's been written for N-terminus with C-terminus
                    auth, entry = writeauths(i, valuesC)
                    authnameC = "authC_" + auth

                    if k == 1:  # before first Nup of a type, C-terminus
                        f.write('\tif (term == "C"):\n')

                    f.write(
                        "\t\t" + authnameC + " = np.array([" + entry + "])\n"
                    )  # coordinates

                    ref = ref + authnameC + ", "

                    # after last Nup of a type, C-terminus
                    if k == count[i[0]]:
                        f.write("\t\tref = np.array([" + ref[:-2] + "])\n\n")
                        k = 0
                        ref = ""
                        mode = "N"  # switch back to C-terminus
                        idlist = []
                    k += 1

            # f.write("\n")

# Cstart = {}
# Nstart = {}
# for i in nups:

#     pattern= "^"+str(i)

#     Cstart[pattern] = valuesC[valuesC.asymID.str.contains(pat=pattern, regex=True)].seqID
#     Nstart[pattern] = valuesN[valuesN.asymID.str.contains(pat=pattern,regex=True)].seqID
