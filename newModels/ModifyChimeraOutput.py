#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 16:43:21 2021

@author: maria
"""

import re
import numpy as np
from warnings import warn

terminus = "N"
rotunit = "multiple"

ChimeraOutput = """
Atom #1.1/G:418@CA 978.520 659.710 677.663
Atom #1.1/M:418@CA 967.247 721.001 720.364
Atom #1.1/S:418@CA 981.755 756.542 750.005
Atom #1.1/Y:418@CA 964.003 697.693 708.018
"""


def CleanCoords(ChimeraOutput, terminus):
    """Used for 5A9Q; rotational unit has to be specified"""

    cleaned = re.sub("#", "", ChimeraOutput)

    # remove whitespaces at start and end of line
    cleaned = re.sub("^ *", "", cleaned, flags=re.MULTILINE)
    cleaned = re.sub(" *$", "", cleaned, flags=re.MULTILINE)

    # remove starting and trailing linebreaks
    while cleaned.index("\n") == 0:
        cleaned = re.sub("^\n", "", cleaned)
        if cleaned.index("\n") != 0:
            cleaned = cleaned[::-1]

    cleaned = re.sub("Atom ", "", cleaned)
    refname = "auth" + terminus + "_"
    cleaned = re.sub("1\../", refname, cleaned)

    cleaned = re.sub(" ", ", ", cleaned)
    cleaned = re.sub(":[0-9]*@CA, ", " = np.array([", cleaned)
    cleaned = re.sub("$", "])", cleaned, flags=re.MULTILINE)

    return cleaned


def CleanCoords2(ChimeraOutput, terminus):
    """Used for 7PEQ; rotational unit does not need to be specified"""
    cleaned = re.sub("", "", ChimeraOutput)

    # remove whitespaces at start and end of line
    cleaned = re.sub("^ *", "", cleaned, flags=re.MULTILINE)
    cleaned = re.sub(" *$", "", cleaned, flags=re.MULTILINE)

    # remove starting and trailing linebreaks
    while cleaned.index("\n") == 0:
        cleaned = re.sub("^\n", "", cleaned)
        if cleaned.index("\n") != 0:
            cleaned = cleaned[::-1]

    cleaned = re.sub("Atom ", "", cleaned)
    refname = "auth" + terminus + "_"
    cleaned = re.sub("/", refname, cleaned)

    cleaned = re.sub(" ", ", ", cleaned)
    cleaned = re.sub(":[0-9]*@CA, ", " = np.array([", cleaned)
    cleaned = re.sub("$", "])", cleaned, flags=re.MULTILINE)

    return cleaned


def summarise_refs(cleaned):
    refslist = re.findall("auth.*=", cleaned)
    refslist = [re.sub("=", "", ref) for ref in refslist]
    refslist = [re.sub(" ", "", ref) for ref in refslist]

    refs = ", ".join([str(elem) for elem in refslist])
    return "ref = np.array([" + refs + "])"


if rotunit == "multiple":
    cleaned = CleanCoords(ChimeraOutput, terminus)
if rotunit == "one":
    cleaned = CleanCoords2(ChimeraOutput, terminus)

refs = summarise_refs(cleaned)
print(cleaned)
print(refs)


def findCentre(nup, nupop):
    """
    nup: Coordinates of a nup of interest
    nupop: Coordinates of the oposite nup of interest.
    For an 8-fold symmetric NPC, the index of an opposite nup is
    + 4.
    e.g. getcrd #1.1/A:20@CA
    and getcrd #1.5/A:20@CA
    """
    # if nup[2] != nupop[2]:
    #     warn("Opposite nups must have the same z coordinates.")
    #     return

    x = np.mean((nup[0], nupop[0]))
    y = np.mean((nup[1], nupop[1]))

    return x, y


nup0 = np.array([1127.765, 714.280, 844.037])
nup1 = np.array([306.035, 712.820, 844.037])
# (716.9000000000001, 713.55)


nup01 = np.array([1035.081, 602.042, 702.602])
nup02 = np.array([398.719, 825.058, 702.602])
# (716.9, 713.55)

nup03 = np.array([828.408, 1031.731, 702.602])
nup07 = np.array([605.392, 395.369, 702.602])
# (716.9000000000001, 713.55)

NupO0 = np.array([1035.081, 602.042, 702.602])
NupO1 = authC_E = np.array([398.719, 825.058, 702.602])

nup7per = np.array([1143.648, 766.967, 1053.077])
nup7per2 = np.array([1045.163, 1421.828, 1052.971])
findCentre(nup7per, nup7per2)


nup107N_00_7r5k = np.array([559.908, 768.369, 1246.970])
nup107N_05_7r5k = np.array([1381.212, 1172.751, 1246.970])
findCentre(nup107N_00_7r5k, nup107N_05_7r5k)

nup107N_00_dilated = np.array([520.991, 919.384, 1252.135])
nup107N_05_dilated = np.array([1466.209, 1067.816, 1252.135])

findCentre(nup107N_00_dilated, nup107N_05_dilated)
