#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:21:05 2021

@author: maria
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import seaborn as sns
from warnings import warn
import matplotlib.animation as animation
import DeformNPC
import Analyse_deformed
from copy import deepcopy
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from scipy.spatial.transform import Rotation as Rot
from sympy.ntheory import primefactors
from scipy.stats import norm
import copy
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from IPython.display import HTML
import NPC

Pos3D = DeformNPC.Pos3D
Sol3D = DeformNPC.Sol3D



class plotOverview:
    def __init__(
        self,
        offsetNPCs,
        NPCs,
        var,
        width=10,
        showforces:bool=False,
        ellipse=False,
        circle=False,
        perRing:bool = False, 
        view:str="front",
        markersizemult=10,
    ):
        """ Plot all NPCs
        :param offsetNPCs: NPC coordinates arranged on a grid 
        :param NPCs: Dictionary containing simulated NPCs and their metadata
        :param var: Dictionary of simulation parameters
        :param width: Width of the figure 
        :param showforces: Show force coordinates, defaults to False
        :param ellipse: Plot fitted ellipses
        :param circles: Plot fitted circles
        :param perRing: Plots circles or ellipses per subcomplex if False, or per Ring if True, defaults to False 
        :param view: Show plot in front- or side view, "front", or "side". Defaults to "front"
        "param markersizemult": Multiplier of plot markersize. Defaults to 10. 
        """
        self.view = view
        self.markersizemult = markersizemult
        if showforces: 
            self.NPCs = NPCs 
            self.var = var
        forcecoords = DeformNPC.Multiple_coord(NPCs["fcoords"], var["symmet"], NPCs["nupIndex"])

        membership = NPCs["z_i_all"] if perRing else NPCs["ringmemall"]

        self.OverviewPlot(
            offsetNPCs,
            forcecoords,
            var["mag"],
            [0],
            membership,
            width,
            showforces=showforces,
            ellipse=ellipse,
            circle=circle,
        )

    def OverviewPlot(
        self,
        offsetNPCs,
        forcecoords,
        mag,
        r,
        membership,
        width,
        showforces=False,
        ellipse=False,
        circle=False,
    ):

        n = int(offsetNPCs[-1, 4] + 1)  # number of NPCs

        if n == 1:
            markersize = self.markersizemult * width
        elif n <= 4:
            markersize = 2 * width
        else:
            markersize = 0.1 * width

        # prepare to colourcode z position
        zs = offsetNPCs[:, 2]
        nup_i = offsetNPCs[:, 3].astype(int)

        zcolour = []
        zcolour.extend(ColourcodeZ(list(zs)))

        # plot
        fig, ax = plt.subplots(1, 1, figsize=(width, width))
        ax.set_title("mag " + str(mag))

        nupmarker = nupMarker()
        nupcolor = nupColor()

        dim2 = 2 if self.view == "side" else 1  # 1 is y-coordinates, 2 is z-coordinates
        for i in np.unique(nup_i):
            x = offsetNPCs[nup_i == i][:, 0]
            y = offsetNPCs[nup_i == i][:, dim2]
            zcolour_subset = [zcolour[j] for j in range(len(zcolour)) if nup_i[j] == i]
            ax.scatter(
                x,
                y,
                c=zcolour_subset,
                s=40 * markersize,
                edgecolors=nupcolor[i],
                marker=nupmarker[i],
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

        if showforces == True:
            offset = NPC.getNPCcoords(self.NPCs, self.var, justoffset = True)
            forceoffset = NPC.getOffsetNPCs(forcecoords, offset = offset)
            ax.scatter(
                forceoffset[:, 0], forceoffset[:, 1], c="blue", s=2 * markersize, marker = "x"
            )

        # fit circle and/or ellipse
        if ellipse:
            NPCs_ellipse = Analyse_deformed.ellipses(offsetNPCs, membership=membership)
            if type(membership[0]) == np.str_:
                for i in range(n):
                    for key in NPCs_ellipse:
                        el = plotEllipse(NPCs_ellipse[key][i])
                        alpha = [float(i) for i in ColourcodeZ(el[:, 2])]
                        ax.scatter(el[:, 0], el[:, 1], lw=2, c="g", alpha=alpha)

            elif type(membership[0]) == np.int64 or type(
                membership[0] == np.float64
            ):  # if one circle fitted per z-ring
                for i in range(len(NPCs_ellipse)):
                    el = plotEllipse(NPCs_ellipse[i])
                    alpha = [float(i) for i in ColourcodeZ(el[:, 2])]
                    ax.scatter(el[:, 0], el[:, 1], lw=1, c="g", alpha=alpha)

        if circle:
            NPCs_circle = Analyse_deformed.circles(offsetNPCs, membership=membership)

            if (
                type(membership[0]) == np.str_
            ):  # if one circle fitted per NPC subcomplexD
                for i in range(n):
                    for key in NPCs_circle:
                        circ = plotCircle(NPCs_circle[key][i])
                        alpha = [float(i) for i in ColourcodeZ(circ[:, 2])]
                        ax.scatter(circ[:, 0], circ[:, 1], lw=2, c="b", alpha=alpha)

            elif type(membership[0]) == np.int64 or type(
                membership[0] == np.float64
            ):  # if one circle fitted per z-ring
                for i in range(len(NPCs_circle)):
                    circ = plotCircle(NPCs_circle[i])

                    alpha = [float(i) for i in ColourcodeZ(circ[:, 2])]
                    if math.isnan(alpha[0]):
                        alpha = 0.7

                    ax.scatter(circ[:, 0], circ[:, 1], lw=1, c="b", alpha=alpha)

        # colourbar for z position
        mincolour = float(min(zcolour))
        maxcolour = float(max(zcolour))

        if mincolour != maxcolour:
            incr = (maxcolour - mincolour) / max(zs)  # increment in colour
            cmap = ListedColormap(
                [
                    str(i)
                    for i in list(np.arange(mincolour, maxcolour + 0.5 * incr, incr))
                ]
            )
            norm = plt.Normalize(min(zs), max(zs))
            fig.colorbar(
                cm.ScalarMappable(norm=norm, cmap=cmap),
                shrink=0.7,
                label="z [nm]",
                ticks=[min(zs), max(zs)],
            )

        ax.set_xticks([])
        ax.set_yticks([])

        barsize = 20
        label = str(barsize) + " nm"
        # label = ""
        scalebar = AnchoredSizeBar(
            ax.transData,
            barsize,
            label,
            "lower right",
            size_vertical=5,
            # bbox_to_anchor = (1050, 70), frameon=False,
            # bbox_to_anchor = (1050, 70),
            frameon=False,
            color=[0.3, 0.3, 0.3],
        )
        ax.add_artist(scalebar)

        ax.axis("scaled")
        fig.tight_layout()


def plotEllipse(NPC):
    return generate_ellipse_by_vectors(NPC["Ce"], NPC["minv"], NPC["majv"])


def plotCircle(NPC):
    (
        r,
        _,
        _,
        _,
        C,
        _,
        normal1,
        normal0,
    ) = NPC  # radius, centre, vectors collinear with circle-plane
    return generate_circle_by_vectors(C, r, normal0, normal1)


def generate_circle_by_vectors(C, r, normal0, normal1):
    t = np.linspace(0, 2 * np.pi, 150)
    P_circle = (
        r * np.cos(t)[:, np.newaxis] * normal0
        + r * np.sin(t)[:, np.newaxis] * normal1
        + C
    )
    return P_circle


def generate_ellipse_by_vectors(C, minv, majv):
    t = np.linspace(0, 2 * np.pi, 100)
    P_ellipse = np.cos(t)[:, np.newaxis] * minv + np.sin(t)[:, np.newaxis] * majv + C
    return P_ellipse


def ColourcodeZ(z, darkest=0.2, brightest=0.7):
    """colourcode z, smaller values are shown darker"""
    return [str(i) for i in np.interp(z, (min(z), max(z)), (darkest, brightest))]


def shifttilt_all_t(pos3D, nFrames, centre, tilt):
    pos3Dc = copy.deepcopy(pos3D)

    if type(tilt) is not type(None):
        for t in range(nFrames):  # tilt for all timeframes
            tiltrot = Rot.from_rotvec(tilt)  #
            pos3Dc[t] = tiltrot.apply(pos3Dc[t])

    for t in range(nFrames):  # shift for all timeframes
        pos3Dc[t, :, 0] += centre[0]
        pos3Dc[t, :, 1] += centre[1]

    return pos3Dc


def plotDetail(
    NPCscoords,
    NPCs:dict,
    var:dict,
    index:int=0,
    width=4,
    mode:str="3D",
    showforces:bool=False,
    trajectory:bool=False,
):
    """ Generate a detail plot of an NPC with given index
    :param NPCscoords: Coordinates of all NPCs
    :param NPCs: Dictionary containing simulated NPCs and their metadata
    :param var: Dictionary of simulation parameters
    :param index: Index of NPC for which the plot should be generated
    :param mode: "3D" or "2D", default: 3D: 2D or 3D plot
    :type mode: str
    :param showforces: Show deforming forces. Will only be correct if tilt and shift have not been applied
    :param trajectories: Show trajectories of nodes over time. Will only be correct if tilt and shift have not been applied
    """

    NoneOrTilt = lambda a, i: a[i] if type(a) is not type(None) else None

    plotNPC = NPCscoords[NPCscoords[:, 4] == index]
    tiltnucv, tiltcytv = DeformNPC.tiltvectors(var["kappa"], var["n"], var["seed"])
    shiftNuc, shiftCyt = DeformNPC.shiftvectors(
        var["shiftsigma"], var["n"], var["seed"]
    )
    tiltnucv_i = NoneOrTilt(
        tiltnucv, index
    )  # tiltnucv[NPCi] if type(tiltnucv) is not type(None) else None
    tiltcytv_i = NoneOrTilt(tiltcytv, index)

    if mode == "2D":
        Plot2D(
            plotNPC,
            NPCs["NPCs"][index],
            NPCs["zexp"],
            NPCs["nupIndex"][index],
            var["symmet"],
            var["nConnect"],
            tiltnucv_i,
            tiltcytv_i,
            shiftNuc[index],
            shiftCyt[index],
            forces=NPCs["fcoords"][index],
            showforces=showforces,
            legend=False,
            trajectory=trajectory,
            colourcode=True,
            springs=False,
            anchorsprings=False,
            width=width,
        )

    elif mode == "3D":
        membership = NPCs["ringmemall"][NPCscoords[:, 4] == index]
        tiltnucv_i = NoneOrTilt(
            tiltnucv, index
        )  # tiltnucv[NPCi] if type(tiltnucv) is not type(None) else None
        tiltcytv_i = NoneOrTilt(tiltcytv, index)
        Plot3D(
            plotNPC,
            NPCs["NPCs"][index],
            NPCs["nupIndex"][index],
            var["nup"],
            var["term"],
            NPCs["isnuc"],
            var["symmet"],
            tiltnucv_i,
            tiltcytv_i,
            shiftNuc[index],
            shiftCyt[index],
            NPCs["fmags"][index],
            NPCs["fcoords"][index],
            membership,
            showforces=showforces,
            trajectory=trajectory,
            viewFrame=-1,
            width=width,
        )


def nupMarker():
    return ["o", "s", "p", "*", "X", "^"]


def nupColor():
    return ["gray", "cyan", "magenta", "green", "purple", "blue"]


def Plot2D(
    plotNPC,
    solution,
    z,
    nup_i,
    symmet,
    nConnect,
    tiltnucv,
    tiltcytv,
    shiftNuc,
    shiftCyt,
    linestyle="-",
    trajectory=True,
    colourcode=True,
    springs=False,
    anchorsprings=False,
    markersize=20,
    forces=0,
    showforces=False,
    legend=False,
    width=5,
):
    """
    solution: Output of solve_ivp
    symmet: number of nodes per ring
    nConnect: number of neighbours connected on each side per node
    linestyle (default: "-"): Linestyle in 1st plot
    legend (default: False): Whether to show a legend in the 1st plot
    colourcode (default: True): colourcodes trajectory with velocity
    colourbar (default: True): Plots colourbar in 2nd plot if True and if colourcode is True
    mainmarkercolor: Colour of nodes in 2nd plot
    """

    nRings = len(z)
    # plt.rcParams.update({'font.size': 25})

    fig, ax = plt.subplots(1, 1, figsize=(width, width))  # 16, 16
    viewFrame = -1  # 0 is the first frame, -1 is the last frame
    mainmarkercolor = ColourcodeZ(z)
    marker = nupMarker()
    edgecolor = nupColor()

    for i in range(nRings):
        centre = shiftNuc if z[i] < np.mean(z) else shiftCyt

        nFrames = len(solution[i].t)  # Nodes at last timestep
        pos3D, vel3D = copy.deepcopy(Sol3D(solution[i]))

        tilt = tiltnucv if z[i] < np.mean(z) else tiltcytv
        pos3D = shifttilt_all_t(pos3D, nFrames, centre, tilt)
        nup = int(nup_i[i])

        ax.plot(
            pos3D[viewFrame, :symmet, 0],
            pos3D[viewFrame, :symmet, 1],
            linestyle="",
            marker=marker[nup],
            color=edgecolor[nup],
            markerfacecolor=mainmarkercolor[i],
            markersize=markersize,
            zorder=50,
            label=str(round(z[i], 1)),
        )

        if anchorsprings:
            ax.plot(centre[0], centre[1], marker="o", color="gray", markersize=15)
            for j in range(symmet):
                ax.plot(
                    (pos3D[viewFrame, j, 0], centre[0]),
                    (pos3D[viewFrame, j, 1], centre[1]),
                    linestyle="-",
                    marker="",
                    color="lightgray",
                    zorder=2,
                    alpha=0.5,
                )

        # circumferential springs
        if springs:
            for ni in range(1, nConnect + 1):  # neighbours to connect to
                for j in range(symmet):  # node to connect from
                    ax.plot(
                        pos3D[viewFrame, (j, (j + ni) % symmet), 0],
                        pos3D[viewFrame, (j, (j + ni) % symmet), 1],
                        linestyle="-",
                        marker="",
                        color="lightgray",
                        zorder=3,
                        alpha=0.5,
                    )  # , linewidth = 5)

        if trajectory:
            if colourcode:  # Colourcoded trajectory
                ### colourcoding velocities
                normvel = np.zeros((nFrames, symmet))  # nFrames, node

                for j in range(symmet):
                    for frame in range(nFrames):
                        normvel[frame, j] = np.linalg.norm(
                            [vel3D[frame, j, 0], vel3D[frame, j, 1], vel3D[frame, j, 2]]
                        )

                norm = plt.Normalize(normvel.min(), normvel.max())

                #####trajectory colorcoded for velocity
                for j in range(symmet):
                    points = pos3D[:, j, :].reshape(-1, 1, 3)  # 3D done
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)

                    lc = LineCollection(
                        segments[:, :, :2],
                        cmap="plasma",
                        norm=norm,
                        zorder=4,
                        linewidth=3.0,
                    )
                    lc.set_array(normvel[:, j])
                    line = ax.add_collection(lc)

            else:  # monochrome trajectory
                for j in range(symmet):
                    ax.plot(pos3D[:, j, 0], pos3D[:, j, 1], color="blue", linestyle="-")

        ### Force vectors
        if showforces and type(forces) != int:
            forces2d = forces[i]
            for j in range(symmet):
                ax.arrow(
                    x=pos3D[0, j, 0],  # frame, node, x or y
                    y=pos3D[0, j, 1],
                    dx=(forces2d[j, 0] - pos3D[0, j, 0]),
                    dy=(forces2d[j, 1] - pos3D[0, j, 1]),
                    width=0.75,
                    color="blue",
                    zorder=101,
                    alpha=1,
                )

    if trajectory and colourcode and legend:
        # fig.legend(bbox_to_anchor=(0.1,-0.025), loc="lower left", ncol = 4, title = "z [nm]")#loc="best")
        axcb = fig.colorbar(line, ax=ax, shrink=0.7, aspect=50, ticks=[0], pad=0.01)
        axcb.set_label("velocity (a.u.)")

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color="0.7",
            label="CR",
            markerfacecolor="0.7",
            markeredgecolor="0.7",
            mew=6,
            markersize=20,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color="0.2",
            label="NR",
            mew=6,
            markersize=20,
        ),
        Line2D(
            [0],
            [0],
            marker="$\\rightarrow$",
            linestyle="",
            color="blue",
            label="force",
            markersize=30,
        ),
        Line2D(
            [0],
            [0],
            marker="",
            linestyle="-",
            color="lightgray",
            label="spring",
            alpha=0.5,
        ),
    ]  # Line2D([0], [0], marker='', linestyle="-", color='m', label='fit', linewidth=2)]

    ax.legend(
        handles=legend_elements,
        loc="upper right",
        bbox_to_anchor=(0.9, -0.1),
        ncol=4,
        fontsize=24,
        frameon=False,
    )

    ax.axis("scaled")
    ax.set(xlabel="x (nm)", ylabel="y (nm)")
    plt.tight_layout()


def positionVStime(NPCs, index=0, width=3, legend=False):
    """
    NPCs: Dictionary containing info of all NPCs
    index: Index of NPC for which plot should be generated
    width: Width of each subplot
    legend: Show ledged. Default False
    """
    symmet = len(
        NPCs["fmags"][index][0]
    )  # length of array of deforming forces on NPC i ring 0. Should be equivalent to symmetry
    XYoverTime(NPCs["NPCs"][index], symmet=symmet, width=width, legend=legend)


def XYoverTime(solution, symmet, width=5, legend=False):
    """x and y positions over time"""

    nRings = len(solution)

    # determin number of rows and colums in final figure. One plot per NPC ring
    l = 2 - nRings % 2
    rows, cols = sorted((int((nRings / l)), l))

    fig, ax = plt.subplots(rows, cols, figsize=(width * rows, width * cols))
    palette = sns.color_palette("hsv", 2 * symmet)

    for ring in range(nRings):
        for i in range(symmet):
            if nRings > 1:
                ax = ax.flatten()
            labelx = labely = None
            if ring == 0:
                labelx = "x" + str(i)
                labely = "y" + str(i)
                labelz = "z" + str(i)

            if nRings > 1:
                ax[ring].plot(
                    solution[ring].t,
                    Pos3D(solution[ring])[:, i, 0],
                    label=labelx,
                    linestyle="-",
                    color=palette[i * 2],
                )
                ax[ring].plot(
                    solution[ring].t,
                    Pos3D(solution[ring])[:, i, 1],
                    label=labely,
                    linestyle="--",
                    color=palette[i * 2],
                )
                ax[ring].plot(
                    solution[ring].t,
                    Pos3D(solution[ring])[:, i, 2],
                    label=labelz,
                    linestyle="-.",
                    color=palette[i * 2],
                )
                ax[ring].set_title("ring " + str(ring))
                ax[ring].set(xlabel="t (a.u.)")
                ax[ring].set(ylabel="change in x or y [nm]")

            else:
                ax.plot(
                    solution[ring].t,
                    Pos3D(solution[ring])[:, i, 0],
                    label=labelx,
                    linestyle="-",
                    color=palette[i * 2],
                )
                ax.plot(
                    solution[ring].t,
                    Pos3D(solution[ring])[:, i, 1],
                    label=labely,
                    linestyle="--",
                    color=palette[i * 2],
                )
                ax.plot(
                    solution[ring].t,
                    Pos3D(solution[ring])[:, i, 2],
                    label=labelz,
                    linestyle="-.",
                    color=palette[i * 2],
                )
                ax.set_title("ring " + str(ring))
                ax.set(xlabel="t (a.u.)")
                ax.set(ylabel="change in x or y [nm]")

    if legend:
        fig.legend(bbox_to_anchor=(1, 0.9), loc="upper left")
    plt.tight_layout()

    plt.show()


def Plot3D(
    plotNPC,
    solution,
    nup_i,
    nupname,
    termname,
    isnuc,
    symmet,
    tiltnucv,
    tiltcytv,
    shiftNuc,
    shiftCyt,
    fmags,
    fcoords,
    membership,
    showforces=False,
    trajectory=False,
    viewFrame=-1,
    colour=["black", "gray"],
    width=5,
):
    """viewFrame: 0 is first frame, -1 is last frame"""
    fig = plt.figure(figsize=(width, width))
    ax = fig.add_subplot(111, projection="3d")

    linewidth = 3
    nRings = len(nup_i)

    marker = nupMarker()
    color = nupColor()

    zs = plotNPC[:, 2]
    colour = ColourcodeZ(zs)

    nupold = None
    for ring in range(nRings):
        # ax.scatter(Pos3D(solution[ring])[viewFrame, : ,0], Pos3D(solution[ring])[viewFrame, :,1], Pos3D(solution[ring])[viewFrame, :,2], s = 300, c = str(colour[ring]), linewidths = linewidth,  marker = "o")

        x = plotNPC[ring * symmet : (ring + 1) * symmet, 0]
        y = plotNPC[ring * symmet : (ring + 1) * symmet, 1]
        z = zs[ring * symmet : (ring + 1) * symmet]

        nup = int(nup_i[ring])
        if ring > 0:
            nupold = int(nup_i[ring - 1])
        label = nupname[nup] + " " + termname[nup] if nup != nupold else ""

        ax.scatter(
            x,
            y,
            z,
            s=300,
            c=colour[ring * symmet : (ring + 1) * symmet],
            edgecolors=color[nup],
            linewidths=linewidth,
            marker=marker[nup],
            label=label,
        )
    ax.legend()

    if showforces:
        for ring in range(nRings):
            x = plotNPC[ring * symmet : (ring + 1) * symmet, 0]
            y = plotNPC[ring * symmet : (ring + 1) * symmet, 1]
            z = zs[ring * symmet : (ring + 1) * symmet]

            for node in range(symmet):
                ax.quiver(
                    x[node],
                    y[node],
                    z[node],
                    fcoords[ring][node][0],
                    fcoords[ring][node][1],
                    0,
                    length=fmags[ring][node],
                    normalize=True,
                    linewidth=linewidth,
                    edgecolor="blue",
                )

    ax.set_xlabel("x [nm]", labelpad=30)
    ax.set_ylabel("y [nm]", labelpad=30)
    ax.set_zlabel("z [nm]", labelpad=40, fontsize=40)

    ax.set_xticks([-25, 25])
    ax.set_yticks([-25, 25])
    ax.set_zticks([0, 50])

    ax.set_xlim3d(-60, 60)
    ax.set_ylim3d(-60, 60)
    ax.set_zlim3d(0, 100)

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Trajectories
    if (
        tiltcytv is not None
        or tiltnucv is not None
        or np.any(shiftCyt)
        or np.any(shiftNuc)
    ):
        trajectory = False
    else:
        trajectory = trajectory

    if trajectory:
        ### colourcoding velocities
        for i in range(len(solution)):  # equals number of rings
            nFrames = len(solution[i].t)  # for NPC 0
            normvel = np.zeros((nFrames, symmet))  # nFrames, node

            pos3D, vel3D = copy.deepcopy(Sol3D(solution[i]))  # for NPC 0

            centre = shiftNuc if isnuc[i] else shiftCyt
            tilt = tiltnucv if isnuc[i] else tiltcytv

            pos3D = shifttilt_all_t(pos3D, nFrames, centre, tilt)

            for j in range(symmet):
                for frame in range(nFrames):
                    normvel[frame, j] = np.linalg.norm(
                        [vel3D[frame, j, 0], vel3D[frame, j, 1], vel3D[frame, j, 2]]
                    )

            norm = plt.Normalize(normvel.min(), normvel.max())

            #####trajectory colorcoded for velocity
            for j in range(symmet):
                points = pos3D[:, j, :].reshape(-1, 1, 3)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = Line3DCollection(
                    segments, cmap="plasma", norm=norm, zorder=4, linewidth=3.0
                )
                lc.set_array(normvel[:, j])
                line = ax.add_collection(lc)

        axcb = fig.colorbar(line, ax=ax, shrink=0.7, aspect=50, ticks=[0], pad=0.01)
        axcb.set_label("velocity (a.u.)")

    # color of edges that aren't axes
    ax.xaxis.pane.set_edgecolor("w")
    ax.yaxis.pane.set_edgecolor("w")
    ax.zaxis.pane.set_edgecolor("w")
    ax.set_box_aspect([1, 1, 1])
    ax.grid(False)

    fitellipse = False
    if fitellipse:
        NPC_ellipse = Analyse_deformed.ellipses(plotNPC, membership=membership)
        if type(membership[0]) == np.str_:
            for key in NPC_ellipse:
                el = plotEllipse(NPC_ellipse[key][0])
                ax.plot(el[:, 0], el[:, 1], el[:, 2], "g-", lw=4)

        elif type(membership[0]) == np.int64 or type(membership[0] == np.float64):
            for i in range(len(NPC_ellipse)):
                el = plotEllipse(NPC_ellipse[i])
                ax.plot(el[:, 0], el[:, 1], el[:, 2], "g-", lw=4)

    fitcircle = False
    if fitcircle:
        NPCs_circle = Analyse_deformed.circles(plotNPC, membership=membership)

        if type(membership[0]) == np.str_:
            for key in NPCs_circle:
                circ = plotCircle(NPCs_circle[key][0])
                ax.plot(circ[:, 0], circ[:, 1], circ[:, 2], "b-", lw=4)

        elif type(membership[0]) == np.int64 or type(membership[0] == np.float64):
            for i in range(len(NPCs_circle)):
                circ = plotCircle(NPCs_circle[i])
                ax.plot(circ[:, 0], circ[:, 1], circ[:, 2], "b-", lw=4)

    plt.show()


# %matplotlib qt


def AnimateDetail(NPCs, var, index:int=0, width=8, directory="", name=None, ext=".gif"):
    """ Animate an individual structure 
    :param NPCs: Dictionary containing simulated NPCs and their metadata
    :param var: Dictionary of simulation parameters
    :param index: Index of structure to be animated, default 0. 0-indexed 
    :param width: Width of animation plot 
    :param directory: directory to save animation in 
    :param name: name of animation. No file will be saved if None. Defaults to None 
    :param ext: File extension, ".gif" or ".mp4", defaults to ".gif"
    """
    AnimatedScatter(
        NPCs["NPCs"][index],
        var["nConnect"],
        var["symmet"],
        NPCs["rexp"],
        NPCs["fmags"][index],
        NPCs["zexp"],
        NPCs["nupIndex"][index],
        width,
        directory,
        name,
        ext,
    )


class AnimatedScatter(object):


    def __init__(
        self,
        solution,
        nConnect,
        symmet,
        r,
        fmags,
        z,
        nup_i,
        width,
        directory,
        name,
        ext,
    ):
        """An animated scatter plot using matplotlib.animations.FuncAnimation."""
        self.solution = solution
        self.nRings = len(solution)
        self.anchors = np.zeros(3)  # anchors
        self.zcolour = ColourcodeZ(z)
        self.nup_i = nup_i
        framenumbers = []

        for i in range(self.nRings):  # check framenumbers are consistent for each ring
            framenumbers.append(len(self.solution[i].t))
        if len(set(framenumbers)) != 1:
            warn(
                "Number of timesteps for all ring must be the same in order to animate deformation."
            )
            return

        nframes = len(self.solution[0].t)

        self.nConnect = nConnect
        self.symmet = symmet
        self.xy = self.xydata()
        self.stream = self.data_stream(self.xy)

        # Setup the figure and axes...
        self.axscale = 1.25 * (np.amax(fmags) + max(r))
        self.fig, self.ax = plt.subplots(figsize=(width, width))
        # plt.rcParams.update({'font.size': 20})

        # Then setup FuncAnimation.

        if ext == "HTML":
            HTML(self.ani.to_html5_video())

        # if name: self.ani.save(name + ".mp4", dpi = 250, fps = 50)
        if name:
            self.ani = animation.FuncAnimation(
                self.fig,
                self.update,
                interval=(5000 / nframes),
                init_func=self.setup_plot,
                blit=True,
                save_count=len(solution[0].t),
            )
            self.ani.save(directory + name + ext, fps=50, dpi=80)

        else:
            self.ani = animation.FuncAnimation(
                self.fig,
                self.update,
                interval=(5000 / nframes),
                init_func=self.setup_plot,
                blit=True,
                save_count=len(solution[0].t),
            )

        plt.show()

    def xydata(self):
        xy = []
        for ring in range(self.nRings):
            xy.append(
                Pos3D(self.solution[ring])[:, np.append(np.arange(self.symmet), 0)]
            )

        return xy  # [xy[0][:,:,:2]] # back to 2D

    def setup_plot(self):
        """Initial drawing of the scatter plot."""

        self.lines = []
        nupcolors = nupColor()
        nupmarker = nupMarker()

        for i in range(
            int(
                self.nRings
                + self.nRings * self.symmet
                + self.symmet * self.nRings * self.nConnect
            )
        ):
            if i < self.nRings:  # all ring-nodes
                self.lobj = self.ax.plot(
                    [],
                    [],
                    marker=nupmarker[self.nup_i[i]],
                    color=self.zcolour[
                        i
                    ],  # str(np.mean(self.xy[i][0][:,2])), # does not accept color
                    markeredgecolor=nupcolors[self.nup_i[i]],
                    linestyle="",
                    markersize=20,
                )
            elif (
                i >= self.nRings and i <= self.nRings * self.symmet
            ):  # ring(s) to anchor
                self.lobj = self.ax.plot(
                    [], [], marker="", color="orange", linestyle="", alpha=0.1, zorder=0
                )  # anchor
            else:  # circumferential springs
                self.lobj = self.ax.plot(
                    [],
                    [],
                    marker="",
                    color="gray",
                    alpha=0.1,
                    linewidth=4,
                    linestyle="-",
                )

            self.lines.append(self.lobj)

        self.ax.axis("scaled")
        self.ax.set(xlabel="x (nm)", ylabel="y (nm)")
        self.ax.axis([-self.axscale, self.axscale, -self.axscale, self.axscale])

        return [
            self.lines[i][0]
            for i in range(
                int(self.nRings * 2 + self.symmet * self.nRings * self.nConnect)
            )
        ]

    def data_stream(self, pos):
        x = np.zeros((self.symmet + 1, self.nRings))
        y = np.zeros((self.symmet + 1, self.nRings))

        while True:
            for i in range(len(self.xy[0])):  # number of timesteps
                for ring in range(self.nRings):
                    x[:, ring] = self.xy[ring][i][
                        :, 0
                    ]  # x coordinates of ring at timestep i
                    y[:, ring] = self.xy[ring][i][:, 1]
                yield x, y

    def update(self, i):
        """Update the plot."""

        x, y = next(self.stream)

        # Alternating anchor coordinates and x or y coordinates
        xa = np.zeros((2 * self.symmet, self.nRings))
        ya = np.zeros((2 * self.symmet, self.nRings))

        for ring in range(self.nRings):  # alternating coordinates of anchor and node
            for i in range(1, 2 * self.symmet, 2):
                # i = 1,3,5,7,9,11,13,15
                xa[i, ring] = x[int((i - 1) / 2), ring]
                ya[i, ring] = y[int((i - 1) / 2), ring]

        xlist = list(x.T) + list(
            xa.T
        )  # 1st 9 entries are nodes 0-7 and 0, following entries are xa
        ylist = list(y.T) + list(ya.T)

        for i in range(self.nRings):
            self.lines[i][0].set_data(xlist[i], ylist[i])  # i is 0, 1, 2, 3

        for i in range(self.nRings, self.nRings * self.symmet + 1):  # to anchor?
            self.lines[i][0].set_data(
                xlist[2 * (i - 1) : 2 * (i - 1) + 2],
                ylist[2 * (i - 1) : 2 * (i - 1) + 2],
            )

        count = (
            self.nRings + self.nRings * self.symmet
        )  # len(xlist) # for one nRIngs == 1, len(xlist) is 2
        for lnum in range(self.nRings):
            for ni in range(1, self.nConnect + 1):  # neighbours to connect to
                for i in range(self.symmet):  # node to connect from
                    self.lines[count][0].set_data(
                        (xlist[lnum][i], xlist[lnum][(i + ni) % self.symmet]),
                        (ylist[lnum][i], ylist[lnum][(i + ni) % self.symmet]),
                    )
                    count += 1
        # return [self.lines[i][0] for i in range(int(self.nRings*2 + self.symmet*self.nRings*self.nConnect))]
        return [
            self.lines[i][0]
            for i in range(
                int(
                    self.nRings
                    + self.nRings * self.symmet
                    + self.symmet * self.nRings * self.nConnect
                )
            )
        ]


if __name__ == "__main__":
    # a = AnimatedScatter(solution, nRings, nConnect, symmet, r, fmags)

    plt.show()


def AnimateOverview(
    NPCs, offsetNPCs, var, width=8, directory="", name=None, ext=".gif"
):
    """An animated scatter plot using matplotlib.animations.FuncAnimation.
    :param NPCs: Dictionary containing simulated NPCs and their metadata
    :param offsetNPCs: NPC coordinates arranged on a grid 
    :param var: Dictionary of simulation parameters
    :param width: width of animation plot 
    :param directory: directory to save animation in, defaults to "". 
    :param name: name of animation. No file will be saved if None. Defaults to None 
    :param ext: File extension, ".gif" or ".mp4", defaults to ".gif"
    """
    AnimateAll(
        NPCs["NPCs"],
        offsetNPCs,
        NPCs["fcoords"],
        var["symmet"],
        NPCs["rexp"],
        width,
        directory,
        name,
        ext,
    )


class AnimateAll(object):
    def __init__(
        self, NPCs, offsetNPCs, forcecoords, symmet, r, width, directory, name, ext
    ):
        self.n = len(NPCs)
        self.tmax = len(NPCs[0][0].t)
        self.t_all = NPCs[0][0].t
        self.maxr = max(r)
        # self.NPCs = NPCs
        self.symmet = symmet
        self.nRings = len(NPCs[0])
        self.nup_i_all = list(
            offsetNPCs[:, 3]
        )  # indicates which nup each point belongs to

        # self.anccoords = deepcopy(anccoords)
        self.forcecoords = deepcopy(forcecoords)

        NPCsCopy = deepcopy(NPCs)

        # Determine the number of rows and columns needed. The last cells on the grid might stay empty
        self.ncols = math.ceil(np.sqrt(self.n))
        self.nrows = math.ceil(self.n / self.ncols)

        self.nrows = 1 if self.n == 1 else max(primefactors(self.n))
        self.ncols = (
            1 if self.n == 1 else int(self.n / max(primefactors(self.n)))
        )  # math.ceil(self.n/self.ncols)

        # i = 0 # will get updated
        # for row in range(self.ncols):
        #     for col in range(self.nrows):
        #         if (i < self.n):

        #             # self.anccoords[i*self.symmet : (i+1)*self.symmet, 0] += col * 8 * self.maxr
        #             # self.anccoords[i*self.symmet : (i+1)*self.symmet, 1] += row * 8 * self.maxr

        #             self.forcecoords[i*self.symmet : (i+1)*self.symmet, 0] += col * 8 * self.maxr
        #             self.forcecoords[i*self.symmet : (i+1)*self.symmet, 1] += row * 8 * self.maxr
        #             i += 1

        self.xy = self.xydata(NPCsCopy)
        self.stream = self.data_stream(self.xy)

        # Setup the figure and axes...
        # self.axscale = 70
        marginm = 1
        self.xmin = round(min(self.xy[0][:, 0])) - marginm * self.maxr
        self.ymin = round(min(self.xy[0][:, 1])) - marginm * self.maxr
        self.xscale = round(max(self.xy[0][:, 0])) + marginm * self.maxr
        self.yscale = round(max(self.xy[0][:, 1])) + marginm * self.maxr

        self.fig, self.ax = plt.subplots(figsize=(width, width))

        self.ani = animation.FuncAnimation(
            self.fig,
            self.update,
            interval=(5000 / self.tmax),
            init_func=self.setup_plot,
            blit=True,
            save_count=self.tmax,
        )

        # if name: self.ani.save(name + ".mp4", dpi = 250, fps = 30)
        if name:
            self.ani.save(directory + name + ext, dpi=80, fps=30)
        if ext == "HTML":
            HTML(self.ani.to_html5_video())
        else:
            plt.show()

    def xydata(self, NPCsCopy):
        xy0 = []

        xy = []
        for npc in range(self.n):
            for ring in range(self.nRings):
                xy0.append(
                    Pos3D(NPCsCopy[npc][ring])
                )  # [:, np.append(np.arange(self.symmet), 0)])

        for frame in range(self.tmax):
            i = 0  # will get updated.
            for col in range(self.ncols):
                for row in range(self.nrows):
                    if i < self.n * self.nRings:
                        for ring in range(self.nRings):
                            xy0[i][frame][:, 0] += col * 3 * self.maxr  # x
                            xy0[i][frame][:, 1] += row * 3 * self.maxr  # y
                            i += 1

        for frame in range(self.tmax):
            xy.append(
                np.reshape(
                    np.array(xy0)[:, frame, :, :],
                    (self.symmet * self.n * self.nRings, 3),
                )
            )

        return xy

    def setup_plot(self):
        # anchors and forces (remains static)

        x, y, t, c = next(self.stream)
        nupcolor = nupColor()
        ec = [
            nupcolor[int(i)] for i in self.nup_i_all
        ]  # edgecolor for points in plot, coded by nup
        # self.ax.scatter(self.anccoords[:,0], self.anccoords[:,1], s = 5, alpha = 0.1, c = "orange")
        # self.ax.scatter(self.forcecoords[:,0], self.forcecoords[:,1], s = 5, alpha = 0.1, c = "cyan")

        # to be updated

        def truncate_colormap(cmap, minval=1, maxval=0.4, n=4):
            new_cmap = colors.LinearSegmentedColormap.from_list(
                "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)),
            )
            return new_cmap

        new_cmap = truncate_colormap(plt.get_cmap("Greys"), n=2 * self.nRings)

        self.scat = self.ax.scatter(
            x,
            y,
            marker="o",
            c=c,
            edgecolor=ec,  # np.random.uniform(0, 1, len(c)),
            s=75,
            label="test",
            cmap=new_cmap,
        )

        self.time_text = self.ax.text(0.02, 0.95, "", transform=self.ax.transAxes)

        self.ax.axis("scaled")
        self.ax.set(xlabel="x (nm)", ylabel="y (nm)")
        self.ax.axis([self.xmin, self.xscale, self.ymin, self.yscale])
        self.time_text.set_text("")
        return self.scat, self.time_text

    def data_stream(self, pos):
        # x = np.zeros(self.symmet * self.n) # why not *nRings?
        # y = np.zeros(self.symmet * self.n)

        while True:
            for frame in range(self.tmax):  # number of timesteps
                #  for ring in range(self.nRings):
                x = self.xy[frame][:, 0]  # x coordinates of ring at timestep i
                y = self.xy[frame][:, 1]
                c = self.colourcodeZ(self.xy[frame][:, 2])
                t = self.t_all[frame]
                # c = list(np.repeat("0.5", len(self.xy[t][:, 0])))#list(str(self.xy[t][:, 2]))
                yield x, y, t, c

    def update(self, i):
        """Update the plot."""

        x, y, t, c = next(self.stream)
        self.scat.set_offsets((np.vstack((x, y)).T))
        self.scat.set_array(np.array(c))
        # roundto = 10
        self.time_text.set_text("t: " + str(int(t)))  # str(roundto * round(t/roundto)))

        return self.scat, self.time_text

    def colourcodeZ(self, z, darkest=0.0, brightest=1):
        """colourcode z, smaller values are shown darker, returns array"""
        return np.interp(z, (min(z), max(z)), (darkest, brightest))


class gethistdata:
    def __init__(
        self, var, NPCs, featuresAll, featuresElAll, featuresel3DAll, width=5, bins=None
    ):
        """
        Plots the distribution of different features for all NPCs. Features are 
        NPC radius, minor/major axis ratio of an ellipse, residual sum of squares of a ring-wise fitted ellipse, 
        difference in tilt in radians between the nucleoplasmic and cytoplasmic rings, determined using ellipses to the respective subcomplexes. 
        shift of nucleoplasmic and cytoplasmic ring, determined using ellipses fitted to the respective subcomplexes. 
        Distance between the centroid of the nucleoplasmic and cytoplasmic ring, determined using ellipses fitted to the respective subcomplexes. 
        :param var: Dictionary of simulation parameters
        :param NPCs: Dictionary containing simulated NPCs and their metadata
        :param featuresAll: features of the NPC determined by fitting circles. From Analyse_deformed.meanfeaturesC()
        :param featuresElAll: features of the NPC determined by fitting ellipses. From Analyse_deformed.meanfeaturesE()
        :param featuresel3DAll: centre and tilt of NPC subcomplexes determined by fitting ellipses. from exportCSV.colfeaturesEl3D. 
        :param width: Width of plot, defaults to 5 
        :param bins: Number of bins to be plotted. If None, the number of bins will be automatically determined. Defaults to None
        """
        self.width = width
        self.bins = bins
        # featuresAll = Analyse_deformed.meanfeaturesC(NPCs, var, circle_allrings)
        # featuresElAll = Analyse_deformed.meanfeaturesE(NPCs, var, ellipse_allrings)
        # _, _, _, _, featuresel3DAll = exportCSV.col_features(NPCs, circle_CRNR, ellipse_CRNR)

        self.rc1 = featuresAll[:, 0]  # circle radius
        self.qe1 = featuresElAll[:, 2]  # minor/major axis ratio
        self.RSSe1 = featuresElAll[:, 4]  # RSS 3D
        self.tilte1 = Analyse_deformed.angles(featuresel3DAll[:, 6:])  # tilt angles
        self.shifte1, self.diste1 = Analyse_deformed.findshiftEl(featuresel3DAll)
        self.featureHist(
            [self.rc1, self.qe1, self.RSSe1, self.tilte1, self.shifte1, self.diste1]
        )

    #   return [rc1, qe1, RSSe1, tilte1, shifte1, diste1]

    def featureHist(self, features1):
        rows = 3
        cols = 2
        name = [
            "radius [nm]",
            "minor/major ellip.",
            "RSS ellip. [nm]",
            "tilt-dif. ellip. [rad]",
            "shift ellip. [nm]",
            "dist ellip [nm]",
        ]

        plt.rcParams.update({"font.size": self.width * 2})

        fig, ax = plt.subplots(
            rows, cols, figsize=(self.width * 0.8, self.width), sharey=False
        )
        # bins = int(len(features1[0])/10) #50

        data = features1[0]

        if self.bins:
            bins = self.bins
        elif data.max() - data.min():
            # determine number of bins via Freedman-Diaconis

            iqr = np.subtract(*np.percentile(data, [75, 25]))
            binwidth = 2 * iqr * (len(data) ** (-1 / 3))
            bins = int((data.max() - data.min()) / binwidth)
            # print(bins)
        else:
            bins = 20

        kwargs = dict(histtype="stepfilled", alpha=0.5, bins=bins, ec="k", linewidth=4)

        for i in range(rows):
            for j in range(cols):
                k = int(i * 2 + j * 1)

                if k == 0:  # set label only for first plot
                    label1 = "$D_{mag}$: 1, $r_\sigma$: 2"
                else:
                    label1 = None

                ax[i, j].hist(features1[k], color="teal", label=label1, **kwargs)

                # ax[i,j].set_ylim([0, 200])
                # ax[i,j].set_yticks([0, 100])
                ax[i, j].set_xlabel(name[k])

        ax[1, 0].set_ylabel("n NPCs")

        def pdffeats(features, i, j, color, label=None, linestyle="-"):
            # i and j are ax indices
            mu, std = norm.fit(features)

            # find area under curve
            values, bns, _ = ax[i, j].hist(features, bins=bins, alpha=0, linewidth=0)
            area = sum(np.diff(bns) * values)

            xmin, xmax = min(bns), max(bns)
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, std) * area  # (xmax-xmin)/bins
            # label = None
            ax[i, j].plot(x, p, color, linewidth=5, linestyle=linestyle, label=label)

        color1 = "darkcyan"
        pdffeats(features1[0], 0, 0, color1)
        pdffeats(features1[1], 0, 1, color1)
        pdffeats(features1[5], 2, 1, color1)

        fig.legend(loc="upper right")
        # fig.legend.get_frame().set_facecolor((0, 0, 1, 0.1))
        # handles, labels = ax.get_legend_handles_labels()
        # fig.legend(handles, labels, loc='upper center')

        plt.tight_layout()
