#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 18:13:18 2021

@author: maria
"""
import math
import numpy as np
import circle_fit as cf
from warnings import warn
from itertools import compress


def ellipses(NPCscoords, membership=None):
    """Fit ellipses to the NPC or its subcomponents 
    
    :param NPCscoords: Coordinates of the NPC 
    :param membership: information to which part of the NPC ellipses should be fitted. 
    if membership = NPCs["ringmemall"], ellipses will be fitted to subcomplexes. 
    if mempership = NPCs["z_i_all"], ellipses will be fitted to rings. 
    if membership = None, ellipses will be fitted to the whole NPC. Defaults to None. 
    :type membership: NoneType, or numpy.ndarray
    :return: features of fitted ellipses per NPC 
    :rtype: list for membership = NPCs["z_i_all"] or membership = None, dict (for membership = NPCs["ringmemall"])
    
    """
    n = len(np.unique(NPCscoords[:, -1]))  # number of NPCs
    iplus1 = int(NPCscoords[-1][-1] + 1)  # index+1 of last or only NPC

    if n != iplus1 and n != 1:
        warn("Indexing of NPCs off.")
        return

    def squaresum(NPC_i, x0, y0, ellipsefeatures_i, z=None, ring=None):
        _, dist = align2ellipse(NPC_i, x0, y0, ellipsefeatures_i, z, ring)
        return sum([j**2 for j in dist])

    def roundfeatures(ellipsefeatures, z=None, ring=None):
        if z != None:
            suffix = "_z" + str(z)
        elif ring != None:
            suffix = "_" + str(ring)
        else:
            suffix = ""

        ellipsefeatures["el_major" + suffix] = round(
            ellipsefeatures["el_major" + suffix], ndigits=1
        )
        ellipsefeatures["el_minor" + suffix] = round(
            ellipsefeatures["el_minor" + suffix], ndigits=1
        )
        ellipsefeatures["el_ssum" + suffix] = round(
            ellipsefeatures["el_ssum" + suffix], ndigits=1
        )
        ellipsefeatures["el_ssumXY" + suffix] = round(
            ellipsefeatures["el_ssumXY" + suffix], ndigits=1
        )
        ellipsefeatures["el_ssumZ" + suffix] = round(
            ellipsefeatures["el_ssumZ" + suffix], ndigits=1
        )
        ellipsefeatures["el_q" + suffix] = round(
            ellipsefeatures["el_q" + suffix], ndigits=2
        )
        ellipsefeatures["el_rot" + suffix] = round(
            ellipsefeatures["el_rot" + suffix], ndigits=3
        )
        return ellipsefeatures

    if str(membership) == "None":
        ellipsefeatures = [None] * n
        squaresumsXY = [None] * n
        squaresumsZ = [None] * n
        squaresums = [None] * n

        i0 = 0
        for i in range(iplus1 - n, iplus1):
            NPC_i = NPCscoords[NPCscoords[:, 4] == i]
            NPC_xy, x0, y0, ellipsefeatures[i0] = fitEllipse3D(NPC_i)
            squaresumsXY[i0] = squaresum(NPC_xy, x0, y0, ellipsefeatures[i0])
            squaresumsZ[i0] = sum(NPC_xy[:, 2] ** 2)
            squaresums[i0] = squaresumsXY[i0] + squaresumsZ[i0]

            ellipsefeatures[i0].update(
                {
                    "el_ssum": squaresums[i0],
                    "el_ssumXY": squaresumsXY[i0],
                    "el_ssumZ": squaresumsZ[i0],
                }
            )

            ellipsefeatures[i0] = roundfeatures(ellipsefeatures[i0])
            i0 += 1

        return ellipsefeatures

    # Fit ellipse to NPC subcomplexes
    elif type(membership[0]) == np.str_:

        def ringmem(name):
            part = NPCscoords[membership == name]
            elPart = [None] * n
            squaresumsXY = [None] * n
            squaresumsZ = [None] * n
            squaresums = [None] * n

            i0 = 0
            for i in range(iplus1 - n, iplus1):
                NPC_part_i = part[part[:, 4] == i]
                NPC_xy, x0, y0, elPart[i0] = fitEllipse3D(NPC_part_i)

                squaresumsXY[i0] = squaresum(NPC_xy, x0, y0, elPart[i0])
                squaresumsZ[i0] = sum(NPC_xy[:, 2] ** 2)
                squaresums[i0] = squaresumsXY[i0] + squaresumsZ[i0]

                elPart[i0].update(
                    {
                        "el_ssum": squaresums[i0],
                        "el_ssumXY": squaresumsXY[i0],
                        "el_ssumZ": squaresumsZ[i0],
                    }
                )

                elPart[i0] = roundfeatures(elPart[i0])

                i0 += 1
            return elPart

        partDictEl = makePartDict(membership, ringmem)

        return partDictEl  # [ellipseCR, ellipseNR]

    elif type(membership[0]) == np.int64 or type(membership[0] == np.float64):
        ellipseAll = [None] * (membership[-1] - membership[0] + 1)

        squaresumsXY = [None] * (membership[-1] - membership[0] + 1)
        squaresumsZ = [None] * (membership[-1] - membership[0] + 1)
        squaresums = [None] * (membership[-1] - membership[0] + 1)
        i0 = 0

        for i in range(membership[0], membership[-1] + 1):
            NPC_xy, x0, y0, ellipseAll[i0] = fitEllipse3D(NPCscoords[membership == i])
            squaresumsXY[i0] = squaresum(NPC_xy, x0, y0, ellipseAll[i0])
            squaresumsZ[i0] = sum(NPC_xy[:, 2] ** 2)
            squaresums[i0] = squaresumsXY[i0] + squaresumsZ[i0]

            ellipseAll[i0].update(
                {
                    "el_ssum": squaresums[i0],
                    "el_ssumXY": squaresumsXY[i0],
                    "el_ssumZ": squaresumsZ[i0],
                }
            )

            ellipseAll[i0] = roundfeatures(ellipseAll[i0])
            i0 += 1

        return ellipseAll


def fitEllipse3D(NPC, z=None, ring=None):
    NPC_centered, V, Ce = fitPlane(
        NPC[:, :3]
    )  # find plane that minimises distance to points

    normal0 = V[0, :]
    normal1 = V[1, :]
    normal2 = V[2, :]

    basisZ = np.array([0, 0, 1])

    # planar rotation only necessary if SVD detects varation in Z
    rotate = False if np.array_equal(basisZ, abs(normal2)) else True
    NPC_xy = rodrigues_rot(NPC_centered, normal2, basisZ) if rotate else NPC_centered  #

    ### fit ellipse parameters to planar NPC
    X = np.array([NPC_xy[:, 0]]).T  # X coordinates NPC
    Y = np.array([NPC_xy[:, 1]]).T  # Y

    # Formulate and solve the least squares problem ||A0x - b ||^2
    A0 = np.hstack([X**2, X * Y, Y**2, X, Y])
    b0 = np.ones_like(X)

    x = np.linalg.lstsq(A0, b0, rcond=None)[0].squeeze()

    A, B, C, D, E = x
    F = -1

    a = -np.sqrt(
        2
        * (A * E**2 + C * D**2 - B * D * E + (B**2 - 4 * A * C) * F)
        * ((A + C) - np.sqrt((A - C) ** 2 + B**2))
    ) / (B**2 - 4 * A * C)

    b = -np.sqrt(
        2
        * (A * E**2 + C * D**2 - B * D * E + (B**2 - 4 * A * C) * F)
        * ((A + C) + np.sqrt((A - C) ** 2 + B**2))
    ) / (B**2 - 4 * A * C)

    x0 = (2 * C * D - B * E) / (
        B**2 - 4 * A * C
    )  # centre of ellipse fitted to NPC_xy. Not necessarily 0, 0
    y0 = (2 * A * E - B * D) / (B**2 - 4 * A * C)

    minor = min(a, b)  # length of minor axis
    major = max(a, b)  # length of major axis

    if B != 0:
        rho = np.arctan(1 / B * (C - A - np.sqrt((A - C) ** 2 + B**2)))
    elif A < C:
        rho = 0
    elif A > C:
        rho = 0.5 * np.pi
    if A < 0:
        rho -= 0.5 * np.pi
        if rho < -0.5 * np.pi:
            rho += np.pi

    q = minor / major

    # points where minor and major axes intercept ellipse fitted to NPC_xy
    xmajor = major * np.cos(rho)
    ymajor = major * np.sin(rho)
    xminor = minor * np.cos(rho - np.pi / 2)
    yminor = minor * np.sin(rho - np.pi / 2)

    # tilted same as NPC, but not offset. Offset applied later via Ce

    if rotate:
        majv = rodrigues_rot(np.array([xmajor, ymajor, 0]), basisZ, normal2)
        minv = rodrigues_rot(np.array([xminor, yminor, 0]), basisZ, normal2)
        Ce = Ce + rodrigues_rot(np.array([x0, y0, 0]), basisZ, normal2)

    else:
        majv = np.array([xmajor, ymajor, 0])
        minv = np.array([xminor, yminor, 0])

    Ce = Ce.flatten()

    keyList = [
        "el_minor",
        "el_major",
        "el_q",
        "el_rot",
        "minv",
        "majv",
        "Ce",
        "normal2",
        "normal1",
        "normal0",
    ]
    valueList = [minor, major, q, rho, minv, majv, Ce, normal2, normal1, normal0]

    if z != None:
        keyList = [i + "_z" + str(z) for i in keyList]
    elif ring != None:
        keyList = [i + "_" + str(ring) for i in keyList]

    zip_iterator = zip(keyList, valueList)
    return NPC_xy, x0, y0, dict(zip_iterator)


def dist2ellipse(semi_major, semi_minor, p):
    # https://github.com/0xfaded/ellipse_demo/issues/1
    # project all points onto 1st quadrant
    px = abs(p[0])
    py = abs(p[1])

    t = math.pi / 4

    a = semi_major
    b = semi_minor

    for _ in range(0, 3):  # 3 iterations
        x = a * math.cos(t)
        y = b * math.sin(t)

        ex = (a * a - b * b) * math.cos(t) ** 3 / a
        ey = (b * b - a * a) * math.sin(t) ** 3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)

        delta_c = r * math.asin((rx * qy - ry * qx) / (r * q))
        delta_t = delta_c / math.sqrt(a * a + b * b - x * x - y * y)

        t += delta_t
        t = min(math.pi / 2, max(0, t))

    closestpoint = (
        math.copysign(a * math.cos(t), p[0]),
        math.copysign(b * math.sin(t), p[1]),
    )
    return closestpoint, np.linalg.norm(p - closestpoint)


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def align2ellipse(NPC, x0, y0, elfeatures, z=None, ring=None):
    # elfeatures = elfeatures[0]
    #    "x0" : x0, "y0" : y0, "minor" : minor, "major" : major, "q": q, "rho" : rho}

    if z != None:
        suffix = "_z" + str(z)
    elif ring != None:
        suffix = "_" + str(ring)
    else:
        suffix = ""

    minor = elfeatures["el_minor" + suffix]
    major = elfeatures["el_major" + suffix]
    rho = elfeatures["el_rot" + suffix]

    nodes = len(NPC)
    NPCrot = np.array(NPC)

    # rotate around centre of fitted ellipse so that the major axis aligns with
    # the y axis
    for i in range(nodes):
        NPCrot[i][0:2] = rotate([x0, y0], NPC[i][0:2], -rho)

    # offset so that the centre is at 0
    NPCrot[:, 0] -= np.mean(NPCrot[:, 0])
    NPCrot[:, 1] -= np.mean(NPCrot[:, 1])

    dist = [None] * nodes
    closestpoints = [None] * nodes

    for i in range(nodes):
        closestpoints[i], dist[i] = dist2ellipse(major, minor, NPCrot[i][0:2])

    return closestpoints, dist


#########################################################################


def circles(NPCscoords, membership=None):
    """Fit circles to the NPC or its subcomponents 
    
    :param NPCscoords: Coordinates of the NPC 
    :param membership: information to which part of the NPC circles should be fitted. 
    if membership = NPCs["ringmemall"], circles will be fitted to subcomplexes. 
    if mempership = NPCs["z_i_all"], circles will be fitted to rings. 
    if membership = None, circles will be fitted to the whole NPC. Defaults to None. 
    :type membership: NoneType, or numpy.ndarray
    :return: features of fitted circles per NPC 
    :rtype: list for membership = NPCs["z_i_all"] or membership = None, dict (for membership = NPCs["ringmemall"])
    
    """
    n = len(np.unique(NPCscoords[:, -1]))  # number of NPCs
    iplus1 = int(NPCscoords[-1][-1] + 1)  # index+1 of last or only NPC

    if n != iplus1 and n != 1:
        warn("Indexing of NPCs off.")
        return

    # if str(ringmode) == "None":
    if str(membership) == "None":
        circlefeatures = [None] * n

        i0 = 0
        for i in range(iplus1 - n, iplus1):
            circlefeatures[i0] = fitCircle3D(
                NPCscoords[NPCscoords[:, 4] == i]
            )  # TODO: 4
            i0 += 1
        return circlefeatures  # {"circlefeatures": circlefeatures} #TODO: Maybe change back

    # if ringmode == "CRNR":
    elif type(membership[0]) == np.str_:
        #
        def ringmem(name):  # name is "CR", "NR", BRCR", "BRNR" or "IR"
            part = NPCscoords[membership == name]
            circlePart = [None] * n
            i0 = 0
            for i in range(iplus1 - n, iplus1):
                circlePart[i0] = fitCircle3D(part[part[:, 4] == i])
                i0 += 1
            return circlePart

        partDict = makePartDict(membership, ringmem)

        return partDict

    # if ringmode == "z":
    elif type(membership[0]) == np.int64 or type(membership[0] == np.float64):
        circlefeaturesZ = [None] * (membership[-1] - membership[0] + 1)
        i0 = 0

        for i in range(membership[0], membership[-1] + 1):
            circlefeaturesZ[i0] = fitCircle3D(NPCscoords[membership == i])
            i0 += 1

        return circlefeaturesZ  # {"circleAll" : circleAll}


def makePartDict(membership, ringmem):
    partDict = {}
    if "CR" in membership:
        partDict["CR"] = ringmem("CR")

    if "NR" in membership:
        partDict["NR"] = ringmem("NR")

    if "BRCR" in membership:
        partDict["BRCR"] = ringmem("BRCR")

    if "BRNR" in membership:
        partDict["BRNR"] = ringmem("BRNR")

    if "IR" in membership:
        partDict["IR"] = ringmem("IR")

    return partDict


def fitCircle3D(NPC):
    NPC = NPC[:, :3]  ##############

    NPC_centered, V, C = fitPlane(NPC)

    normal0 = V[0, :]
    normal1 = V[1, :]
    normal2 = V[2, :]  # if no variation in z, normal2= [0, 0, 1]

    basisZ = np.array([0, 0, 1])
    rotate = False if np.array_equal(basisZ, normal2) else True

    NPC_xy = rodrigues_rot(NPC_centered, normal2, basisZ) if rotate else NPC_centered  #

    zsqsum = sum(NPC_xy[:, 2] ** 2)  # square sum of distances from circle-plane in z

    xc, yc, r, residual = cf.least_squares_circle(
        NPC_xy[:, :2]
    )  # xcentre, ycentre, radius

    sqsum = zsqsum + residual

    if rotate:
        C = C + rodrigues_rot(
            np.array([xc, yc, 0]), basisZ, normal2
        )  # Transform circle center back to 3D coords

    C = C.flatten()

    return [
        r,
        sqsum,
        residual,
        zsqsum,
        C,
        normal2,
        normal1,
        normal0,
    ]  # [xc, yc, np.round(r, 2), residual] # TODO: change back


def fitPlane(P):
    P_mean = P.mean(axis=0)  # centre
    P_centered = P - P_mean  # 0-centre
    U, s, V = np.linalg.svd(P_centered)
    return P_centered, V, P_mean


def rodrigues_rot(P, n0, n1):
    # If P is only 1d array (coords of single point), fix it to be matrix
    if P.ndim == 1:
        P = P[np.newaxis, :]

    # Get vector of rotation k and angle theta
    n0 = n0 / np.linalg.norm(n0)
    n1 = n1 / np.linalg.norm(n1)

    k = np.cross(n0, n1)
    k = k / np.linalg.norm(k)
    theta = np.arccos(np.dot(n0, n1))

    # Compute rotated points
    P_rot = np.zeros((len(P), 3))
    for i in range(len(P)):
        P_rot[i] = (
            P[i] * np.cos(theta)
            + np.cross(k, P[i]) * np.sin(theta)
            + k * np.dot(k, P[i]) * (1 - np.cos(theta))
        )

    return P_rot


def fitCircle(NPC):
    xc, yc, r, residual = cf.least_squares_circle(
        NPC
    )  # xcentre, ycentre, radius #TODO: cf only works in 2D
    xc = round(xc, ndigits=1)
    yc = round(yc, ndigits=1)
    # r = round(r, ndigits=1) TODO
    residual = round(residual, ndigits=1)
    return [
        xc,
        yc,
        r,
        residual,
    ]  # [xc, yc, np.round(r, 2), residual] # TODO: change back


def meanfeaturesC(NPCs, var, circle_allrings:list, full:bool = False):
    """
    :param NPCs: Dictionary containing simulated NPCs and their metadata
    :param var: Dictionary of simulation parameters
    :param circle_allrings: features of circles fitted to each NPC ring 
    :param full: indicates whether to also return centreAll and tiltAll, default False 
    :return: featuresAll, shape (n NPCs, 4)
    featuresAll[:,0]: mean NPC circle radius 
    featuresAll[:,1]: mean sum of squared errors (SSE) of circles fitted per NPC 
    featuresAll[:, 2]: mean SSE in lateral direction of circles fitted per NPC 
    featuresAll[:, 3]: mean SSE in axial direction of circles fitted per NPC 
    optional: 
    centreAll: shape (n NPCs, 3): mean centre per NPC, determined by fitting circles 
    tiltAll:  shape (n NPCs, 3): mean tilt of NPC, determined by fitting circles 
    """
    z_i_ring = np.repeat(
        [i for i in range(var["n"])], len(NPCs["zexp"])
    )  # index that assigns each ring to an NPC
    n = len(np.unique(z_i_ring))  # n NPCs
    featuresAll = np.zeros((n, 4))  # 4 because radius, SSE, SSE lateral, SSE axial



    for i in range(n):
        circle_NPC_n = list(
            compress(circle_allrings, z_i_ring == i)
        )  # circle_allrings, but only values for NPC with index i
        range_n_rings_in_i = range(
            len(circle_NPC_n)
        )  # range of numbers of rings within NPC i
        featuresAll[i] = np.array(
            [
                np.mean([circle_NPC_n[j][k] for j in range_n_rings_in_i])
                for k in range(4)
            ]
        )  # r, sqsum, residual, zsqsum
        
        
    if full: 
        centreAll = np.zeros((n, 3))  # 3 because 3D
        tiltAll = np.zeros((n, 3))  # 3 because 3D     
        for i in range(n):
                
            
            centreAll[i] = np.mean(
                [circle_NPC_n[j][4] for j in range_n_rings_in_i], axis=0
            )  # centre per NPC
            tiltAll[i] = np.mean(
                [circle_NPC_n[j][5] for j in range_n_rings_in_i], axis=0
            )  # tilt per NPC
        return featuresAll, centreAll, tiltAll
    
    return featuresAll

def meanfeaturesE(NPCs, var, ellipse_allrings:list, el_name=False, full = False):
    """
    :param NPCs: Dictionary containing simulated NPCs and their metadata
    :param var: Dictionary of simulation parameters
    :param ellipse_allrings: features of ellipses fitted to each NPC ring
    :type ellipse_allrings: list of dictionaries
    :param el_name: list of feature names, corrsponding to keys in ellipse_allring. 
    ["el_major", "el_minor", "el_q" "el_rot", "el_ssum", "el_ssumXY", "el_ssumZ"] if False. Defaults to False. 
    :return: 
    featuresElAll: Mean features per NPC determined by averaging ring-wise fitted ellipses. 
        shape (nNPCs, n features), n features is 7 by default, corresponding to el_name: 
        length major axis, length minor axis, ratio minor/major, azimuthal angle ellipse, 
        sum of square errors, sum of square errors in lateral direction, sum of square errors in axial direction
    centreElAll: shape (n NPCs, 3): mean centre per NPC, determined by fitting ellipses 
    tiltElAll:  shape (n NPCs, 3): mean tilt of NPC, determined by fitting ellipses 
    """
    z_i_ring = np.repeat([i for i in range(var["n"])], len(NPCs["zexp"]))

    if el_name == False:
        el_name = [
            "el_major",
            "el_minor",
            "el_q",
            "el_rot",
            "el_ssum",
            "el_ssumXY",
            "el_ssumZ",
        ]
    n = len(np.unique(z_i_ring))  # n NPCs
    featuresElAll = np.zeros((n, len(el_name)))


    for i in range(n):
        ellipse_NPC_n = list(compress(ellipse_allrings, z_i_ring == i))
        range_n_rings_in_i = range(len(ellipse_NPC_n))

        featuresElAll[i] = np.array(
            [
                np.mean([ellipse_NPC_n[j][i] for j in range_n_rings_in_i])
                for i in el_name
            ]
        )
    
    if full: 
        centreElAll = np.zeros((n, 3))  # 3 because 3D
        tiltElAll = np.zeros((n, 3))  # 3 because 3D    
        for i in range(n):
            centreElAll[i] = np.mean(
                [ellipse_NPC_n[i]["Ce"] for i in range_n_rings_in_i], axis=0
            )
            tiltElAll[i] = np.mean(
                [ellipse_NPC_n[i]["normal2"] for i in range_n_rings_in_i], axis=0
            )
    
        return featuresElAll, centreElAll, tiltElAll
    return featuresElAll


def angles(entriestiltAll):
    """Input: n * 6 array of tilt angles, where the first 3 columns are tilt-angles of
    the first sub-complex, and the next 3 colums are tilt-angles of the second sub-complex.
    n is the number of NPCs"""

    def findangle(v1, v2):  # vector 1, vector 2
        uv1 = v1 / np.linalg.norm(v1)
        uv2 = v2 / np.linalg.norm(v2)
        return np.arccos(np.dot(uv1, uv2))

    tiltangles = np.zeros(len(entriestiltAll))

    for i in range(len(entriestiltAll)):
        tiltv1 = entriestiltAll[i][:3]
        tiltv2 = entriestiltAll[i][3:]
        tiltangles[i] = findangle(tiltv1, tiltv2)
        if tiltangles[i] >= np.pi / 2:
            tiltangles[i] = np.pi - tiltangles[i]

    return tiltangles


def findshiftEl(featuresel3DAll):
    c1 = featuresel3DAll[:, :2]
    c2 = featuresel3DAll[:, 3:5]
    d1 = featuresel3DAll[:, 2]
    d2 = featuresel3DAll[:, 5]
    shift = np.linalg.norm(c1 - c2, axis=1)
    dist = d1 - d2
    return shift, dist
