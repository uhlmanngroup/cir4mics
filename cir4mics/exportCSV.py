#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:21:00 2023

@author: maria
"""

from datetime import date
import uuid
import csv
import math
import Analyse_deformed
import numpy as np
from itertools import compress


def names(var):
    def makeNames(mean, sigma, meanName, sigmaName, forTitle=True):
        meanStr = str(meanName) + "_" + str(mean) + "_" if mean else ""
        sigmaStr = str(sigmaName) + "_" + str(sigma) + "_" if sigma else ""
        return meanStr + sigmaStr

    today = date.today()
    rStr = makeNames(var["rnew"], var["rsigma"], "r_mean", "r_sigma")
    dStr = makeNames(var["dnew"], var["dsigma"], "dist_mean", "dist_sigma")
    thetaStr = makeNames(
        var["thetanew"], var["thetasigma"], "twist_ang_mean", "twist_ang_sigma"
    )
    elongStr = makeNames(
        var["elliptnew"], var["elliptsigma"], "elong_mean", "elong_sigma"
    )
    symmetStr = "_symmet_" + str(var["symmet"]) + "_" if var["symmet"] != 8 else ""
    kappaStr = makeNames(var["kappa"], None, "kappa", None)
    shiftStr = makeNames(var["shiftsigma"], None, "ring-shift_sigma", None)

    nup_index = ["Ch_" + str(i) for i in list(range(len(var["term"])))]
    thetaref = (
        "theta_ref:_" + var["nup"][0] + "_" + var["term"][0] + "_" if var["rel"] else ""
    )

    channels = [
        sub[item]
        for item in range(len(var["term"]))
        for sub in [nup_index, list(var["nup"]), list(var["term"])]
    ]
    channels = ["_".join(channels)][0]
    modelname = "model_" + str(var["model"])

    # name = (channels+ "_x" + str(n) + modelname
    #     + "_deformmag_" + str(mag) + "_" +  symmetStr + rStr + dStr + thetaStr + thetaref +
    #     elongStr + kappaStr + shiftStr + "seed_" + str(seed) + "_" + today.strftime("%d-%b-%Y"))

    uuid_name = str(uuid.uuid1())
    name = (
        channels
        + "_x"
        + str(var["n"])
        + "_"
        + modelname
        + "_deformmagxy_"
        + str(var["mag"])
        + "_"
        + uuid_name
        + "_"
        + today.strftime("%d-%b-%Y")
    )

    return {
        "rStr": rStr,
        "dStr": dStr,
        "thetaStr": thetaStr,
        "thetaref": thetaref,
        "elongStr": elongStr,
        "symmetStr": symmetStr,
        "kappaStr": kappaStr,
        "shiftStr": shiftStr,
        "channels": channels,
        "modelname": modelname,
    }, name


class MakeCSV:
    def __init__(self, var, NPCs, printNPC, nameDict, name, data_dir):
        """Export CSV file. printNPC is either offsetNPCs or NPCscoords.
        Choose NPCscoords if photophysics simulations software introduces offsets to NPCs, i.e. arranges them on a grid
        """

        self.csvpath = data_dir + name + ".csv"
        with open(self.csvpath, "w", newline="") as csvfile:
            writer1 = csv.writer(
                csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )
            writer1.writerow(["x", "y", "z", "channel", "particle"])

            for i in range(len(printNPC)):
                # writer1.writerow([printNPC[i,0], printNPC[i,1], printNPC[i,2], printNPC[i,3].astype(int)])
                writer1.writerow(
                    list(printNPC[i, :3])
                    + [printNPC[i, 3].astype(int)]
                    + [printNPC[i, 4].astype(int)]
                )

        # corresponding metadata file
        f = open(data_dir + name + "_metadata.txt", "x")

        # f.write(str(term) + "_" + str(nup) + "_" + str(model) + "\n")
        f.write(nameDict["channels"] + "\n")
        f.write(nameDict["modelname"] + "\n")

        f.write("random seed: " + str(var["seed"]) + "\n")
        f.write("\nmodel layout: ###########\n")

        # f.write("kr = 0.7 \n")
        f.write("nConnect: " + str(var["nConnect"]) + "\n")

        f.write("\nNPC parameters: ##################")

        f.write("\nIrregular deformation parameters: ###########\n")

        f.write("deformation magnitude, xy: " + str(var["mag"]) + "\n")
        f.write("deformation magnitude, z: " + str(var["zmag"]) + "\n")
        f.write("sigma multiplier " + str(var["sigmamult"]) + "\n")

        f.write("\nGeometric NPC properties (input values): ###########\n")
        f.write("symmetry: " + str(var["symmet"]) + "\n")
        f.write(
            "Nup index per ring: (corresponding to channel): "
            + str(NPCs["nupIndex"][0])
            + "\n"
        )
        f.write(
            "Ring membership to subcomplex: " + str(list(NPCs["ringmember"])) + "\n"
        )
        f.write("expected radii per ring: " + str(NPCs["rexp"]) + "\n")
        f.write("expected z per ring: " + str(NPCs["zexp"]) + "\n")
        f.write(
            "expected ring angles per ring (rad): " + str(NPCs["ringAnglesExp"]) + "\n"
        )
        f.write("twist angle input (rad): " + str(NPCs["thetaold"]) + "\n")
        f.write(
            "Nearest rotational unit, NR - CR (degrees): "
            + str(math.degrees(NPCs["theta_offset"]))
            + "\n"
        )
        f.write("\n")

        if nameDict["rStr"]:
            f.write("changed radii: " + nameDict["rStr"].replace("_", " ") + "\n")
        if nameDict["dStr"]:
            f.write("changed distance: " + nameDict["dStr"].replace("_", " ") + "\n")
        if nameDict["elongStr"]:
            f.write(
                "minor/major axis input: "
                + nameDict["elongStr"].replace("_", " ")
                + "\n"
            )
        if nameDict["thetaref"]:
            f.write(
                "ref nup for ring-angles and rotation, "
                + nameDict["thetaref"].replace("_", " ")
                + "\n"
            )
        if nameDict["thetaStr"]:
            f.write(
                "twist angles (rad): " + nameDict["thetaStr"].replace("_", " ") + "\n"
            )
        if nameDict["kappaStr"]:
            f.write(
                "kappa of von mises-fisher distribution: " + str(var["kappa"]) + "\n"
            )
        if nameDict["shiftStr"]:
            f.write(
                "shift of rings (sigma of 0-centered normal, nm): "
                + str(var["shiftsigma"])
                + "\n"
            )
        f.close()


class featuresCSV:
    def __init__(self, NPCs, var, name, circle_allrings, ellipse_allrings, data_dir):
        c_name = [
            "c_r",
            "c_sse",
            "c_sselat",
            "c_sseax",
            "c_x0",
            "c_y0",
            "c_z0",
            "c_tilt_x",
            "c_tilt_y",
            "c_tilt_z",
        ]
        el_name = [
            "el_major",
            "el_minor",
            "el_q",
            "el_rot",
            "el_ssum",
            "el_ssumXY",
            "el_ssumZ",
        ]

        self.csvpath = data_dir + "NPC_features_" + name + ".csv"
        with open(self.csvpath, "w", newline="") as csvfile:
            writer2 = csv.writer(
                csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )

            writer2.writerow(
                ["NPC"]
                + c_name
                + el_name
                + ["ce_x0", "ce_y0", "cez0"]
                + ["e_tilt_x", "e_tilt_y", "e_tilt_z"]
            )

            featuresAll, centreAll, tiltAll = Analyse_deformed.meanfeaturesC(
                NPCs, var, circle_allrings, full = True
            )
            featuresElAll, centreElAll, tiltElAll = Analyse_deformed.meanfeaturesE(
                NPCs, var, ellipse_allrings, el_name, full = True
            )

            for i in range(var["n"]):
                writer2.writerow(
                    [i]
                    + list(featuresAll[i])
                    + list(centreAll[i])
                    + list(tiltAll[i])
                    + list(featuresElAll[i])
                    + list(centreElAll[i])
                    + list(tiltElAll[i])
                )


def col_features(NPCs, circle_CRNR, ellipse_CRNR):
    """features per subcomplex"""
    _, idx = np.unique(NPCs["ringmember"], return_index=True)
    memberof = list(
        NPCs["ringmember"][np.sort(idx)]
    )  # list of subcomplexes, e.g. ['CR', 'NR']
    dim = 3  # dimensions

    feat1D = ["r", "SSE", "SSE_l", "SSE_a"]  # 1D features circle
    featEl = [
        "el_minor",
        "el_major",
        "el_q",
        "el_rot",
        "el_ssum",
        "el_ssumXY",
        "el_ssumZ",
        "Ce",
        "normal2",
    ]
    nfEl = len(featEl)

    subs = len(memberof)
    n = len(NPCs["NPCs"])
    entriesAll = np.zeros((n, len(feat1D) * subs))
    entriescenAll = np.zeros((n, dim * subs))
    entriestiltAll = np.zeros((n, dim * subs))
    featureselAll = np.zeros((n, (nfEl - 2) * subs))
    featuresel3DAll = np.zeros((n, 2 * dim * subs))

    for NPC in range(n):
        entriesAll[NPC] = np.array(
            [circle_CRNR[sub][NPC][i] for i in range(len(feat1D)) for sub in memberof]
        )  # e.g. CR, NR, CR, NR ...
        entriescenAll[NPC] = np.array(
            [
                circle_CRNR[sub][NPC][4][i]
                for sub in memberof
                for i in range(len(circle_CRNR[sub][NPC][4]))
            ]
        )
        entriestiltAll[NPC] = np.array(
            [
                circle_CRNR[sub][NPC][5][i]
                for sub in memberof
                for i in range(len(circle_CRNR[sub][NPC][5]))
            ]
        )

        entriesel = [
            ellipse_CRNR[sub][NPC][featEl[i]] for i in range(nfEl) for sub in memberof
        ]  # entries alterate between subcomplexes
        featureselAll[NPC] = np.array(
            entriesel[: len(memberof) * 7]
        )  # 'el_minor', 'el_major', 'el_q', 'el_rot', 'el_ssum', 'el_ssumXY', 'el_ssumZ'
        featuresel3DAll[NPC] = np.array(
            entriesel[len(memberof) * 7 :]
        ).flatten()  # 'Ce', 'normal2', alternating by subcomplex

    return entriesAll, entriescenAll, entriestiltAll, featureselAll, featuresel3DAll


class featuresCSV_subcomplex:
    def __init__(self, NPCs, circle_CRNR, ellipse_CRNR, name, data_dir):
        _, idx = np.unique(NPCs["ringmember"], return_index=True)
        memberof = list(NPCs["ringmember"][np.sort(idx)])

        col = ["r", "SSE", "SSE_l", "SSE_a"]
        col2 = ["Cx", "Cy", "Cz"]
        col3 = ["Tx", "Ty", "Tz"]

        columns = [c + "_" + mem for c in col for mem in memberof]
        cen = [c + "_" + mem for mem in memberof for c in col2]
        tilt = [c + "_" + mem for mem in memberof for c in col3]

        colel = [
            "el_minor",
            "el_major",
            "el_q",
            "el_rot",
            "el_ssum",
            "el_ssumXY",
            "el_ssumZ",
            "Ce",
            "normal2",
        ]
        columnsel = [c + "_" + mem for c in colel[:7] for mem in memberof]
        cenel = ["el_" + c + "_" + mem for mem in memberof for c in col2]
        tiltel = ["el_ " + c + "_" + mem for mem in memberof for c in col3]

        self.csvpath = data_dir + "NPC_features_subcomplex_" + name + ".csv"
        with open(self.csvpath, "w", newline="") as csvfile:  # TODO: change directory
            writer2 = csv.writer(
                csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )

            writer2.writerow(
                ["NPC"] + columns + cen + tilt + columnsel + cenel + tiltel
            )

            (
                entriesAll,
                entriescenAll,
                entriestiltAll,
                featureselAll,
                featuresel3DAll,
            ) = col_features(NPCs, circle_CRNR, ellipse_CRNR)

            n = len(NPCs["NPCs"])

            for NPC in range(n):
                writer2.writerow(
                    [NPC]
                    + list(entriesAll[NPC])
                    + list(entriescenAll[NPC])
                    + list(entriestiltAll[NPC])
                    + list(featureselAll[NPC])
                    + list(featuresel3DAll[NPC])
                )


class featuresCSV_rings:
    def __init__(self, NPCs, var, name, data_dir, circle_allrings, ellipse_allrings):
        flattenvalues = lambda NPC_n, key, outputrange: np.array(
            [NPC_n[j][key] for j in outputrange]
        ).flatten()

        nfeat = 6  # 6 because: r, sse, sse_l, sse_a, centre, tilt
        nfeatel = 9
        col = ["r", "SSE", "SSE_l", "SSE_a"]
        col2 = ["Cx", "Cy", "Cz"]
        col3 = ["Tx", "Ty", "Tz"]

        NP = [i.replace("up", "") for i in var["nup"]]  # shorten nup names
        z = NPCs["zexp"]
        nupIndex = NPCs["nupIndex"]
        colnames = [
            cl + "_" + NP[nupIndex[0][j]] + "_" + str(round(z[j])) + "nm"
            for cl in col
            for j in range(len(z))
        ]

        cen = [
            cl + "_" + NP[nupIndex[0][j]] + "_" + str(round(z[j])) + "nm"
            for j in range(len(z))
            for cl in col2
        ]
        tilt = [
            cl + "_" + NP[nupIndex[0][j]] + "_" + str(round(z[j])) + "nm"
            for j in range(len(z))
            for cl in col3
        ]
        colnames2 = cen + tilt

        col_e = [
            "el_minor",
            "el_major",
            "el_q",
            "el_rot",
            "el_ssum",
            "el_ssumXY",
            "el_ssumZ",
            "Ce",
            "normal2",
        ]

        colel = [
            "el_minor",
            "el_major",
            "el_q",
            "el_rot",
            "el_ssum",
            "el_ssumXY",
            "el_ssumZ",
        ]
        colel2 = ["eCx", "eCy", "eCz"]
        colel3 = ["eTx", "eTy", "eTz"]

        z_i_ring = np.repeat([i for i in range(var["n"])], len(NPCs["zexp"]))

        colnamesel = [
            cl + "_" + NP[nupIndex[0][j]] + "_" + str(round(z[j])) + "nm"
            for cl in colel
            for j in range(len(z))
        ]
        cenel = [
            cl + "_" + NP[nupIndex[0][j]] + "_" + str(round(z[j])) + "nm"
            for j in range(len(z))
            for cl in colel2
        ]
        tiltel = [
            cl + "_" + NP[nupIndex[0][j]] + "_" + str(round(z[j])) + "nm"
            for j in range(len(z))
            for cl in colel3
        ]
        colnamesel2 = cenel + tiltel

        self.csvpath = data_dir + "NPC_features_rings_" + name + ".csv"
        with open(self.csvpath, "w", newline="") as csvfile:  # TODO: change directory
            writer2 = csv.writer(
                csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )

            writer2.writerow(["NPC"] + colnames + colnames2 + colnamesel + colnamesel2)

            for i in range(var["n"]):
                circle_NPC_n = list(compress(circle_allrings, z_i_ring == i))
                ellipse_NPC_n = list(compress(ellipse_allrings, z_i_ring == i))

                range_n_rings_in_i = range(len(circle_NPC_n))

                circleRings = []
                ellipseRings = []
                for k in range(nfeat):
                    circleRings.append(
                        flattenvalues(circle_NPC_n, k, range_n_rings_in_i)
                    )

                for k in range(nfeatel):
                    ellipseRings.append(
                        flattenvalues(ellipse_NPC_n, col_e[k], range_n_rings_in_i)
                    )

                features = np.array(
                    circleRings[: nfeat - 2]
                ).flatten()  # -2 because last two entries are 3D vectors (tilt and centre)
                features3D = np.array(circleRings[nfeat - 2 :]).flatten()

                featuresel = np.array(ellipseRings[: nfeatel - 2]).flatten()
                featuresel3D = np.array(ellipseRings[nfeatel - 2 :]).flatten()

                writer2.writerow(
                    [i]
                    + list(features)
                    + list(features3D)
                    + list(featuresel)
                    + list(featuresel3D)
                )


def previewCSV(namespace, lines=5):
    count = 0
    with open(namespace) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            print(row)
            count += 1
            if count > lines:
                break
