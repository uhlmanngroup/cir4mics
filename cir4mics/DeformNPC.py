#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:36:41 2021

@author: maria
"""

import numpy as np
import numpy.matlib
from scipy.integrate import solve_ivp
from scipy.linalg import circulant, null_space
from warnings import warn
from sklearn.gaussian_process.kernels import RBF
import random
import sys
from Nups_Info import SelectNup
from numba import njit
import Analyse_deformed
import math

np.seterr(
    divide="ignore", invalid="ignore"
)  # To stop printing "RuntimeWarning: invalid value encountered in divide"

# import line_profiler
# profile = line_profiler.LineProfiler()


class DeformNPC:
    def __init__(
        self,
        r,
        ringAngles,
        z,
        symmet=8,
        elliptical=False,
        mag=0,
        zmag=0,
        sigmamult=0.5,
        nConnect=2,
        damp=1,
        kr=0.7,
        tlast=20,
        step=0.25,
        kmult=1,
        seed=None,
    ):
        """
        Models deformed NPCs using solve_ivp based on a simple, rotationally symmetric node and spring model with deforming forces applied in xy direcion.
        ### Input ###
        ## Geometric parameters ##
        nRings: Number of NPC rings
        r: Radius of NPC rings assuming 8-fold symmetry. Must be a list of length nRings
        ringAngles: Rotational angle offset of NPC rings. Must be a list of length nRings
        z: z position of NPC rings. Stays unchanged. Must be a list of length nRings
        symmet: Symmetry of the NPC (Default 8 )

        ## Influencing forces/offsets ##
        elliptical: Minor/major axis ratio of an ellipse. Generates elongated NPCs by sampling elongating forces. Default False. Not equivalent to minor/major axis ratio of output NPC
        mag: Magnitude of deformation in x-y  [nm]. Number represents 3 standard deviations -> 99.7 % of forces on a node lie within this range
        zmag: Magnitude of offset (not using springs) in z
        Sigmamult: influences RBF kernel from which forces are sampled. Multiple of NPC radius. Default 0.5
        nConnect: Number of connected neighbour nodes in clock-wise and anti-clockwise direction
        damp: Damping of springs. By default 1
        kr: spring constant of anchor springs. Default 0.7. Stabilises against drift

        ## Time steps ##
        tlast: last computed time-step, default 20. Can be higher for more complex models/deformations
        step: Time-steps for which output is saved. Default for ever 0.25 timesteps. Higher value results in smoother output but does not change end-result

        ## Randomisation ##
        Seed: Number between 1 and inf. Random seed for reporducibility.

        ### Relevant output ###
        solution: solve_ivp output for all NPC rings
        initcoords: starting coordinates
        fmag: list of magnitude of forces for all NPC rings
        r: Radius of NPC rings corrected for symmetry (i.e. increasing with symmetry)
        z: z position of NPC rings

        """

        self.seed = seed
        self.solution = []
        self.symmet = symmet  # Rotational symmetry of the NPC
        self.initcoords = []  # starting coordinates of nodes
        self.z = z  # z coordinates (stays unchanged)
        self.elliptical = elliptical
        self.geodesic = np.zeros(len(ringAngles))
        self.nRings = len(self.z)

        self.tlast = tlast
        tspan = [0, self.tlast]
        self.teval = np.arange(0, self.tlast, step)  # needed for animation function
        # teval = None

        Lrests = []  # circulant matrix of resting spring lengths
        Ks = (
            []
        )  # circulant matrix of spring constants # TODO: Way to lower and increase spring constants of entire model?
        y0s = []  # initial coordinates and velocities per node

        self.geodesic = np.asarray(ringAngles) * np.mean(
            r
        )  # calculate geodesic distance of corners (for 8-fold NPC)
        self.ringAngles_corrected = np.zeros(len(ringAngles))

        # correct radius
        mean_r_corrected = self.adjustRadiusStatic(np.mean(r), self.symmet)
        diff_r = r - np.mean(r)
        r_corrected = diff_r + mean_r_corrected

        self.ringAngles_corrected = self.geodesic / np.mean(r_corrected)

        for i in range(self.nRings):
            initcoord = self.initialcoords(
                r_corrected[i], self.z[i], ringAngle=self.ringAngles_corrected[i]
            )
            self.initcoords.append(initcoord)
            Lrest = self.springlengths(initcoord)
            Lrests.append(Lrest)
            Ks.append(self.springconstants(Lrest, kmult=kmult))
            y0s.append(np.concatenate((initcoord.flatten(), np.zeros(3 * self.symmet))))

        self.r = r_corrected
        sigma = np.min(self.r) * sigmamult  # sigma: free parameter for RBF kernel

        self.fmag = self.forcesMultivariateNorm(
            self.initcoords, r_corrected, mag, sigma=sigma
        )  # generate random forces
        self.zoffset = self.forcesMultivariateNorm(
            self.initcoords, r_corrected, zmag, sigma=sigma, addim=2
        )

        # Solve ODE, ring 1 - 4
        self.fcoord = np.zeros(np.shape(self.initcoords))

        anchor = np.zeros(3)
        for i in range(self.nRings):
            anchor[2] = self.z[i]
            self.solution.append(
                solve_ivp(
                    self.npc,
                    tspan,
                    y0s[i],
                    t_eval=self.teval,
                    rtol=1e-1,
                    atol=1e-2,
                    args=(
                        self.initcoords[i],
                        r_corrected[i],
                        Lrests[i],
                        Ks[i],
                        kr,
                        self.fmag[i],
                        damp,
                        nConnect,
                        anchor,
                    ),
                )
            )

            self.fcoord[i] = self.initialcoords(
                r_corrected[i], self.z[i], self.ringAngles_corrected[i], self.fmag[i]
            )

    ### Methods ###
    def npc(self, t, y, initcoords, r, Lrest, K, kr, fmag, damp, nConnect, anchor):
        """
        t: time points
        y: values of the solution at t
        r: radius of NPC ring
        Lrest: Circulant matrix of resting lengths of all springs
        K: Circulant matrix of all radial spring constants
        kr: Spring constants of anchor springs
        fmag: array of forces (length = symmet) to be applied in radial direction to each node
        damp: Damping factor
        nConnect: Number of connected neighbours in cw and ccw direction for each node

        output: solutions at t. x and y components of positions and velocities of each node for each timestep
        """

        return self.npcIntAnc(
            self.symmet, y, anchor, r, fmag, nConnect, K, Lrest, kr, damp
        )

    @staticmethod
    @njit
    def npcIntAnc(symmet, y, anchor, r, fmag, nConnect, K, Lrest, kr, damp):
        v = np.reshape(y[3 * symmet :], (symmet, 3))
        x = np.reshape(y[: 3 * symmet], (symmet, 3))

        F = np.zeros((symmet, 3))  # Forces
        for i in range(symmet):
            F[i, :2] = fmag[i] * x[i, :2] / np.linalg.norm(x[i, :2] - anchor[:2])

        allaccarray = np.zeros(
            (symmet, 3)
        )  # array for accelerations of node 0 - symmet-1

        for i in range(symmet):  # i indicates the reference node
            accarray = np.array(
                [0.0, 0.0, 0.0]
            )  # initiate acceleration array for each node i

            for j in [
                k for k in range(-nConnect, nConnect + 1) if k != 0
            ]:  # j is neighbour nodes -nConnect to +nConnect relative to i, skipping 0 (0=i)
                jnew = (i + j) % symmet
                accarray += (
                    K[i][jnew]
                    * (x[i] - x[jnew])
                    / np.linalg.norm(x[i] - x[jnew])
                    * (Lrest[i][jnew] - np.linalg.norm(x[i] - x[jnew]))
                )

            accarray += (
                kr
                * (x[i] - anchor)
                / np.linalg.norm(x[i] - anchor)
                * (r - np.linalg.norm(x[i] - anchor))
            )  # anchor
            accarray = F[i] + accarray - damp * v[i]  # external force and damping
            allaccarray[i] = accarray

        dxdt = np.concatenate((v.flatten(), allaccarray.flatten()))
        return dxdt

    @staticmethod
    @njit
    def adjustRadiusStatic(
        r8, symmet
    ):  # TODO: Why is this part of the deformation code; could be elsewhere
        """Adjusts radius r with symmetry. No adjustment is made when symmetry is 8. Radius is viewed
        as the length of the symmetric side of an isoceles triangle whose tip (angle alpha) points towards the
        center of the NPC and whose base is the section between two neighbouring nodes at the circumference. Think slice of cake.
        # Input:
        r8: radius of a default 8-fold symmetrical NPC

        ## Output:
        radius of NPC with symmetry equal to symmet (rnew = r8 if symmet = 8)

        """
        alpha = (
            2 * np.pi / symmet
        )  # Angle at the tip of triangular slice (pointing to center of NPC)
        theta = 0.5 * (
            np.pi - alpha
        )  # Either angle at the base of (isosceles) triangular slice
        halfbase = r8 * np.sin(
            np.pi / 8
        )  # half the distance between two corners of an NPC ring (= half the base of triangular slice)
        return halfbase / np.cos(theta)  # new radius

    @staticmethod
    @njit
    def pol2cart(rho, phi):
        """Transforms polar coordinates of a point (rho: radius, phi: angle) to 2D cartesian coordinates."""
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return (x, y)

    @staticmethod
    @njit
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return (rho, phi)

    def initialcoords(self, r, z, ringAngle=0, forces=0):
        """
        Generates cartesian coordinates of the NPC given radius and self.symmet
        ## Input ##
        r: NPC Radius
        ringAngle: Angular offset of ring (default 0)
        forces (optional): input 1D forcevector to view coordinates of forces; used for visualisation

        ## Return values ##
        2D Cartesian coordinates
        """

        if type(forces) != int:
            if len(forces) != self.symmet:
                warn("forces must be 0 or an array with len(self.symmet")

        initcoord = self.initialcoordsStatic(r, z, self.symmet, ringAngle, forces)

        return initcoord

    @staticmethod
    @njit
    def initialcoordsStatic(r, z, symmet, ringAngle=0, forces=0):
        forces = forces * np.ones(
            symmet
        )  # forces is 0 or np.array with len(self.symmet)
        rotAngle = np.linspace(
            0, (symmet - 1) * 2 * np.pi / symmet, symmet
        )  # angles of nodes around NPC
        initcoord = np.zeros((symmet, 3))

        initcoord[:, 0] = (r + forces) * np.cos(rotAngle + ringAngle)
        initcoord[:, 1] = (r + forces) * np.sin(rotAngle + ringAngle)
        initcoord[:, 2] = z

        return initcoord

    def springlengths(self, initcoord):
        """Compute lengths of springs from coordinates and returns circulant matrix"""
        l = np.linalg.norm((initcoord[0, :] - initcoord), axis=1)
        return circulant(l)

    def springconstants(self, Lrest, kmult=1):
        "Returns circulant matrix of spring constants"

        Lscale = Lrest / Lrest[0][1]  # scaled so that shortest edge is 1

        k = np.ones(int(np.floor(self.symmet / 2)))
        k *= kmult

        if self.symmet % 2 == 0:  # if symmet is even
            k[-1] = (
                k[-1] / 2
            )  # springs that connect opposite corners will be double-counted. Spring constant thus halved
            K = circulant(np.append(0, np.append(k, np.flip(k[:-1]))))
        else:  # if symmet is odd
            K = circulant(np.append(0, [k, np.flip(k)]))

        return K / Lscale

    def forcesMultivariateNorm(self, initcoords, r, mag, sigma, addim=0):
        """
        Returns array of Forces that are later applied in radial direction to the NPC corners
        ## Input ##
        *coordring: Initial coordinates of nodes for an arbitrary number of rings.
        mag: Total magnitude of distortion. Default is 50.
        r: radius of all NPC rings
        addim: Additional dimmensions. 0 for radial forces
        ## Returns ##
        For each ring, an array of forces applied to each node
        """

        nodesTotal = self.symmet * self.nRings  # total number of nodes over all rings
        allcoords3D = np.asarray(initcoords).reshape(
            nodesTotal, 3
        )  # separated by NPC rings

        # compute offset for ellipticity if minor/major axis isn't None or 1, and only in 2D
        computeFe = (
            True if self.elliptical and self.elliptical != 1 and addim == 0 else False
        )

        if computeFe:
            Fe = self.ellipse(allcoords3D[:, :2], self.elliptical)

        if mag == 0:
            F = np.zeros(nodesTotal)
        else:
            # project all coordinates onto one circle so that
            # varying radii keep their shape relative to each other better
            phi = np.arctan2(allcoords3D[:, 1], allcoords3D[:, 0])
            allcoords3D[:, 0] = np.mean(r) * np.cos(phi)
            allcoords3D[:, 1] = np.mean(r) * np.sin(phi)

            cov = np.zeros((nodesTotal, nodesTotal))

            kernel = 1.0 * RBF(sigma)
            cov = kernel.__call__(X=allcoords3D)

            var = (mag / 3) ** 2  # 3*SD to Var TODO
            cov = var * cov
            cov = (
                cov + np.identity(nodesTotal) * 1e-6
            )  # jitter matrix to surely make it positive definite

            rng = np.random.default_rng(seed=self.seed)
            u = rng.standard_normal((1 + addim) * nodesTotal)
            L = np.linalg.cholesky(cov)
            F = L @ u[addim * nodesTotal :]

        # if NPCs elliptical and anchors internal:
        if computeFe:
            F = F + Fe

        return np.split(F, self.nRings)

    def ellipse(self, allcoords, elliptical):
        allcoords = allcoords[:, :2]
        p = [i * 2 * np.pi for i in self.r]  # perimeter of NPC circle

        q = elliptical  # ratio of a to b

        a = [
            ((np.sqrt(2) * i) / (2 * np.pi)) / np.sqrt(1 + q**2) for i in p
        ]  # length of axis a given perimeter of the ellipse is rougly p
        b = [q * i for i in a]

        nodes = len(allcoords)
        Fe = np.zeros(np.shape(allcoords))  # forces to deform NPC into ellipse

        polarcoords = np.zeros(np.shape(allcoords))
        n = np.zeros(
            nodes
        )  # scaling factors n[i] that transform allcoords into ellipse coords
        ellipsecoords = np.zeros(np.shape(allcoords))

        # rotate NPC ring randomly between -180 and +180 deg
        for i in range(nodes):
            polarcoords[i] = self.cart2pol(allcoords[i, 0], allcoords[i, 1])

        rng = np.random.default_rng(seed=self.seed)
        rotateby = rng.uniform(-np.pi, np.pi)

        polarcoords[:, 1] += rotateby

        for i in range(nodes):
            ring = int(np.floor(i / self.symmet))
            allcoords[i] = self.pol2cart(
                polarcoords[i, 0], polarcoords[i, 1]
            )  # rotated NPC ring back to cartesian coords
            n[i] = np.sqrt(
                (a[ring] ** 2 * b[ring] ** 2)
                / (
                    b[ring] ** 2 * allcoords[i, 0] ** 2
                    + a[ring] ** 2 * allcoords[i, 1] ** 2
                )
            )  # scaling factors n[i] that transform allcoords into ellipse coords
            ellipsecoords[i] = allcoords[i] * n[i]

        difference = ellipsecoords - allcoords

        sign = np.sign(
            difference[:, 0] * allcoords[:, 0] + difference[:, 1] * allcoords[:, 1]
        )

        Fe = sign * np.linalg.norm(difference, axis=1)

        return Fe


# Multiple NPCs


# DeformNPC()
def MultipleNPC(
    nup,
    term,
    model,
    n=1,
    rel=False,
    rnew=None,
    rsigma=None,
    thetanew=None,
    thetasigma=None,
    dnew=None,
    dsigma=None,
    symmet=8,
    elliptnew=None,
    elliptsigma=None,
    mag=0,
    zmag=0,
    sigmamult=0.5,
    nConnect=2,
    damp=1,
    kr=0.7,
    tlast=20,
    step=0.25,
    seed=None,
    kmult=1,
    **kwargs
):
    """
    Simulate n deformed NPCs using solve_ivp based on a simple, rotationally symmetric node and spring model
    with deforming forces applied in xy direcion and axial offset zmag.
    ### Input ###
    ## Select NPC protein to tag
    nup: One or more (as tuple) nucleoporines to tag
    term: protein terminus to tag, "N", or "C"
    model: PDB model to extract coordinates from: "5A9Q", "7PEQ", "5IJN", "5IJO", "7PER", "7R5K", or "7R5J".
    relative: If True, Defines first listed Nup as reference in multi-channel simulations.

    ## Geometric parameters ##
    rvar: Defines changed mean radius and standard deviation of radii. Computes radii from input-model if not defined
    thetavar: Defines changed mean angle and standard deviation of twist angle between nucleoplasmic and cytoplasmic side. Computes values from input-model if not defined
    dvar: Defines changed mean distance and standard deviation of distance between nucleoplasmic and cytoplasmic side. Computes value from input-model if not defined.
    symmet: Symmetry of the NPC (Default 8 )

    ## Influencing forces/offsets ##
    elliptvar: Mean and standard deviation of minor/major axis ratio of an ellipse. Generates elongated NPCs by sampling elongating forces. Default False. Not equivalent to minor/major axis ratio of output NPC
    mag: Magnitude of deformation in x-y  [nm]. Number represents 3 standard deviations -> 99.7 % of forces on a node lie within this range
    zmag: Magnitude of offset (not using springs) in z
    sigmamult: influences RBF kernel from which forces are sampled. Multiple of NPC radius. Default 0.5
    nConnect: Number of connected neighbour nodes in clock-wise and anti-clockwise direction
    damp: Damping of springs. By default 1
    kr: spring constant of anchor springs. Default 0.7. Stabilises against drift

    ## Time steps ##
    tlast: last computed time-step, default 20. Can be higher for more complex models/deformations
    step: Time-steps for which output is saved. Default for ever 0.25 timesteps. Higher value results in smoother output but does not change end-result

    ## Randomisation ##
    Seed: Number between 1 and inf. Random seed for reporducibility.
    """

    if nConnect > symmet / 2:
        nConnect = int(np.floor(symmet / 2))
        warn(
            "Selected number of neighbours nConnect too large. nConnect has been changed to "
            + str(nConnect)
            + "."
        )

    selectedNup = SelectNup(nup, term, model)

    r = selectedNup.r
    z = selectedNup.z
    ringAngles = selectedNup.ringAngles
    nup_i = selectedNup.nupindex
    ringMember = selectedNup.ringmember
    nRings = len(z)

    if len(r) != nRings or len(ringAngles) != nRings:
        warn("r, ringAngles, and z must be of length nRings: " + str(nRings))

    thetaold = selectedNup.rotAng  # NR - CR
    thetaoffset = selectedNup.rotOffset  # 0, -45 or +45. Only works for 8-fold symmetry

    rold = np.array(r)
    zold = np.array(z)
    # NucSideBool should be all True for only one ring
    isnuc = zold <= np.mean(zold)
    NucSideBool = np.tile(
        np.repeat(isnuc, symmet), n
    )  # node on nuclear side or not (all NPCs)
    ringMembAll = np.tile(np.repeat(ringMember, symmet), n)
    z_i_all = np.repeat(np.arange(n * nRings), symmet)

    if (
        rel
    ):  # define first listed Nup as reference for ringAngles, theta, radius, and distance
        refNup = SelectNup(nup[0], term[0], model)
        ringAngles -= np.mean(ringAngles[nup_i == 0])
        thetaold = refNup.rotAng
        thetaoffset = refNup.rotOffset  # 0, -45 or +45. Only works for 8-fold symmetry
        if rnew != None:
            rnew *= np.mean(rold) / np.mean(refNup.r)

        if dnew != None:
            zref = refNup.z
            midplaneold = np.mean(zold)
            midplaneref = np.mean(zref)
            dabs = np.mean(zold[zold > midplaneold]) - np.mean(zold[zold < midplaneold])
            drel = np.mean(zref[zref > midplaneref]) - np.mean(zref[zref < midplaneref])
            dnew += dabs - drel

    ringAnglesOld = np.array(ringAngles)

    # saved per NPC:
    NPCs = []  # ODE output
    fcoords = []  # absolute coordinates of forces
    fmags = []  # magnitudes of forces
    zoffsets = []  # random offset in z (for variability in z)
    newthetas = []  # new rotational offset
    ellipticals = []  # updated ratio minor/major axis
    nup_is = []

    if seed:
        random.seed(seed)
        seeds = random.sample(range(min(2**32, sys.maxsize)), n)  # , n)
    else:
        seeds = np.full(n, None)

    # expected values: pretends no standard deviation is applied, to export as metadata
    rexp = Change_radius(rold, rnew)  # rsigma, etc, False if not specified
    zexp = Change_dist(zold, dnew)

    if math.isnan(
        thetaold
    ):  # thetaold will be nan if all nups lie on the same z-plane.
        ringAnglesExp = ringAnglesOld
        newminTheta = None  # updated twist angle. No twist angle if all nups lie on the same z-plane
    else:
        ringAnglesExp, _ = Change_rotang(
            ringAnglesOld, thetaold, thetaoffset, zold, thetanew
        )

    # Modification 1 here
    # dnew = 30
    # rnew = 0.8 * np.mean(r)
    # elliptnew = 0.835
    # thetanew = 0.244346
    # #kappanew = 10
    # #shiftnew = 5

    # basis = np.array([0, 8, 16, 24, 32, 40])
    # mags = np.repeat(np.array([0, 1, 5, 10, 15, 20]),8)
    ######

    for i in range(n):  # for each NPC
        # update r, z, angles, and ellipticity
        nup_is.append(nup_i)
        r = Change_radius(rold, rnew, rsigma, seeds[i])
        if not math.isnan(thetaold):  # If Nups lie on different z-planes
            z = Change_dist(zold, dnew, dsigma, seeds[i])
            ringAngles, newminTheta = Change_rotang(
                ringAnglesOld,
                thetaold,
                thetaoffset,
                zold,
                thetanew,
                thetasigma,
                seeds[i],
            )
        elliptical = Change_ellipt(elliptnew, elliptsigma, seeds[i])

        # Modification 2 here
        # if i in basis + 6:
        #    z = Change_dist(zold, dnew, dvar["dsigma"], seeds[i])
        # elif i in basis + 5:
        #    r = Change_radius(rold,  rnew, rvar["rsigma"], seeds[i])
        # elif i in basis + 3:
        #    ringAngles, newminTheta = Change_rotang(ringAnglesOld, thetaold, thetaoffset, zold, thetanew, thetavar["thetasigma"], seeds[i])
        # elif i in basis + 2:
        #    elliptical = Change_ellipt(elliptnew, elliptvar["elliptsigma"], seeds[i])

        # mag = mags[i]
        # zmag = mag/2
        ######

        # deform NPC
        # (self r, ringAngles, z, symmet = 8, elliptical = False,
        #  mag = 0, zmag = 0, sigmamult = 0.5, nConnect = 2, damp = 1,
        #  kr = 0.7, tlast = 20, step = 0.25, seed = None):

        deformNPC_temp = DeformNPC(
            r,
            ringAngles,
            z,
            symmet=symmet,
            elliptical=elliptical,
            mag=mag,
            zmag=zmag,
            sigmamult=sigmamult,
            nConnect=nConnect,
            damp=damp,
            kr=kr,
            tlast=tlast,
            step=step,
            kmult=kmult,
            seed=seeds[i],
        )

        # save output
        NPCs.append(deformNPC_temp.solution)
        fcoords.append(deformNPC_temp.fcoord)  # radial
        fmags.append(deformNPC_temp.fmag)
        zoffsets.append(deformNPC_temp.zoffset)
        newthetas.append(newminTheta)
        ellipticals.append(elliptical)

    # updated expected value for r
    rexp = DeformNPC(
        rexp,
        ringAnglesExp,
        zexp,
        symmet=symmet,
        elliptical=elliptical,
        mag=0,
        zmag=0,
        sigmamult=sigmamult,
        nConnect=nConnect,
        seed=seeds[i],
    ).r

    ringAngles_corrected = (
        deformNPC_temp.ringAngles_corrected
    )  # for NPCs that aren't 8-fold symmetric

    multipleNPCs = {
        "NPCs": NPCs,
        "nupIndex": nup_is,
        "sigmamult": sigmamult,
        "zexp": zexp,
        "fcoords": fcoords,
        "rexp": rexp,
        "ringAnglesExp": ringAnglesExp,
        "ringAngles_corrected": ringAngles_corrected,
        "elliptical": ellipticals,
        "fmags": fmags,
        "nConnect": nConnect,
        "thetaold": thetaold,
        "newthetas": newthetas,
        "theta_offset": thetaoffset,
        "zoffsets": zoffsets,
        "nuclear_side_boolean": NucSideBool,
        "isnuc": isnuc,
        "ringmemall": ringMembAll,
        "ringmember": ringMember,
        "z_i_all": z_i_all,
    }

    return multipleNPCs


def Change_radius(r, rnew=False, rsigma=False, seed=None):
    """
    r: old radii of all rings, np.array
    rnew: new mean radius or False
    sigma: standard deviation of Gaussion of which to sample new mean radii from.
    Gaussian is centred around rnew, if rnew is provided, and around np.mean(r) otherwise
    """

    if str(rnew) == "0" or str(rnew) == "0.0":
        sys.exit("Radius can't be 0")

    rmean = rnew if rnew else np.mean(r)

    if rsigma:
        np.random.seed(seed)
        rmean = np.random.normal(rmean, rsigma)

    return (rmean / np.mean(r)) * r


def Change_dist(zold, dnew=False, dsigma=False, seed=None):
    z = np.array(zold)

    midplane = np.mean(z)
    dist = np.mean(z[z > midplane]) - np.mean(z[z <= midplane])

    if not (str(dnew) == "0" or str(dnew) == "0.0"):
        dnew = dnew if dnew else dist

    if dsigma:
        np.random.seed(seed)
        dnew = np.random.normal(dnew, dsigma)

    ddif = dist - dnew

    z[z > midplane] = z[z > midplane] - ddif

    return z


def Change_rotang(
    ringAnglesOld,
    thetaold,
    thetaoffset,
    zold,
    thetanew=False,
    thetasigma=False,
    seed=None,
):
    midplane = np.mean(zold)
    ringAngles = np.array(ringAnglesOld)

    if not (str(thetanew) == "0" or str(thetanew) == "0.0"):
        thetanew = thetanew if thetanew else thetaold

    if thetasigma:
        np.random.seed(seed)
        thetanew = np.random.normal(thetanew, thetasigma)

    thetadif = thetaold - thetanew
    ringAngles[zold > midplane] += thetadif / 2
    ringAngles[zold <= midplane] -= thetadif / 2

    NRmean = np.mean(ringAngles[zold <= midplane])
    CRmean = np.mean(ringAngles[zold > midplane])

    newminTheta = NRmean - CRmean - thetaoffset

    return ringAngles, newminTheta


def Change_ellipt(elliptnew=False, elliptsigma=False, seed=None):
    elliptnew = elliptnew if elliptnew else 1.0
    if elliptsigma:
        np.random.seed(seed)
        elliptnew = np.random.normal(elliptnew, elliptsigma)
    return elliptnew


def Sol3D(solution):  # TODO: 3D
    """
    input: DeformNPC.solution for a given NPC ring
    output: 2D arrays of position and of velocity of nodes in a ring over time [timeframe, node, dimension (x or y)]
    """
    nFrames = len(solution.t)
    symmet = int(
        len(solution.y) / 6
    )  # /6 because of 3 dim times both velocity and position
    pos1D = solution.y.T[:, : 3 * symmet]  # positions over time
    vel1D = solution.y.T[:, 3 * symmet :]  # velocities over time
    pos3D = np.reshape(pos1D, (nFrames, symmet, 3))
    vel3D = np.reshape(vel1D, (nFrames, symmet, 3))
    return pos3D, vel3D


def Pos3D(solution):  # TODO: 3D
    pos3D, vel3D = Sol3D(solution)
    return pos3D


def MultipleNPCs_coord(
    NPCs,
    zoffsets,
    symmet,
    NucSideBool,
    nupIndex,
    tiltnucv,
    tiltcytv,
    shiftNuc,
    shiftCyt,
    seed=None,
    frame=-1,
):  
    "Input is OdeResults for each NPC in a list. Output is just the final coordinates of each NPC"

    nNPCs = len(NPCs)
    nRings = len(NPCs[0])  
    NPCscoord = np.zeros(
        (nNPCs * nRings * symmet, 5)
    )  # number of nodes, dimensions + label #TODO: more general
    # 0 is first frame, -1 is last frame
    i = 0

    for NPC in range(nNPCs):
        for ring in range(nRings):
            NPCscoord[i * symmet : (i + 1) * symmet] = np.c_[
                (
                    Pos3D(NPCs[NPC][ring])[frame, :],
                    np.repeat(nupIndex[NPC][ring], symmet),
                    np.repeat(NPC, symmet),
                )
            ]
            NPCscoord[i * symmet : (i + 1) * symmet, 2] += zoffsets[NPC][
                ring
            ]  # offset in z
            i += 1  # counts rings for all NPCs

    if not isinstance(tiltnucv, type(None)):  # TODO: tilt and shift do not commute
        for NPC in range(nNPCs):
            NPCscoord = tilt(NPC, NPCscoord, NucSideBool, tiltnucv, tiltcytv)

    # shift rings
    if not isinstance(shiftNuc, type(None)):
        for NPC in range(nNPCs):
            NPCscoord = shift(NPC, NPCscoord, NucSideBool, shiftNuc, shiftCyt)

    return NPCscoord


def shiftvectors(shiftsigma, nNPCs, seed):
    dim = 2  # dimensions for shift vector. 2 and shift only in x-y
    if shiftsigma != None:
        rng = np.random.default_rng(seed)
        shiftNuc = rng.normal(0, shiftsigma, (nNPCs, dim))
        shiftCyt = rng.normal(0, shiftsigma, (nNPCs, dim))
    else:
        shiftNuc = np.zeros([nNPCs, dim])
        shiftCyt = np.zeros([nNPCs, dim])
    return shiftNuc, shiftCyt


def tilt(NPC, NPCscoord, NucSideBool, tiltnucv, tiltcytv):
    """NPC: index of NPC
    NPCscoord: coordinates of all NPCs, not shifted,
    NucSideBool: array, True for each node on the nucleoplasmic side
    tiltnucv: for each NPC, rotational vector defining tilt of nucleoplasmic rings
    tiltcytv: for each NPC, rotational vector defining tilt of cytoplasmic rings
    """
    Nuc_of_NPCBool = np.logical_and(NPCscoord[:, 4] == NPC, NucSideBool)
    Cyt_of_NPCBool = np.logical_and(NPCscoord[:, 4] == NPC, np.invert(NucSideBool))

    NPCscoordsN = np.array(NPCscoord[Nuc_of_NPCBool, :3])
    meanzN = np.mean(
        NPCscoord[Nuc_of_NPCBool, 2]
    )  # mean Z position of NPC nucleoplasmic side
    NPCscoordsN[:, 2] -= meanzN

    NPCscoordsC = np.array(NPCscoord[Cyt_of_NPCBool, :3])
    meanzC = np.mean(NPCscoord[Cyt_of_NPCBool, 2])
    NPCscoordsC[:, 2] -= meanzC

    baseZ = np.array([0, 0, 1])

    NPCscoord[Nuc_of_NPCBool, :3] = Analyse_deformed.rodrigues_rot(
        NPCscoordsN, baseZ, tiltnucv[NPC]
    )
    NPCscoord[Cyt_of_NPCBool, :3] = Analyse_deformed.rodrigues_rot(
        NPCscoordsC, baseZ, tiltcytv[NPC]
    )

    NPCscoord[Nuc_of_NPCBool, 2] += meanzN
    NPCscoord[Cyt_of_NPCBool, 2] += meanzC
    return NPCscoord


def shift(NPC, NPCscoord, NucSideBool, shiftNuc, shiftCyt):
    """NPC: index of NPC
    NPCscoord: coordinates of all NPCs, not shifted,
    NucSideBool: array, True for each node on the nucleoplasmic side
    shiftNuc: for each NPC, vector defining shift of nucleoplasmic rings
    shiftCyt: for each NPC, vector defining shift of cytoplasmic rings
    """

    Nuc_of_NPCBool = np.logical_and(NPCscoord[:, 4] == NPC, NucSideBool)
    Cyt_of_NPCBool = np.logical_and(NPCscoord[:, 4] == NPC, np.invert(NucSideBool))

    NPCscoord[Nuc_of_NPCBool, :2] += shiftNuc[NPC]
    NPCscoord[Cyt_of_NPCBool, :2] += shiftCyt[NPC]
    return NPCscoord


def tiltvectors(kappa, nNPCs, seed):  # TODO: cite
    """
    https://dlwhittenbury.github.io/ds-2-sampling-and-visualising-the-von-mises-fisher-distribution-in-p-dimensions.html
    """
    if kappa == None:
        return None, None
    rng = np.random.default_rng(seed)
    mu = [0, 0, 1]  # standard basis vector in z

    def rand_uniform_hypersphere(N, p):
        v = rng.normal(0, 1, (N, p))
        v = np.divide(v, np.linalg.norm(v, axis=1, keepdims=True))
        return v

    def rand_t_marginal(kappa, p, N):
        # Start of algorithm
        b = (p - 1.0) / (2.0 * kappa + np.sqrt(4.0 * kappa**2 + (p - 1.0) ** 2))
        x0 = (1.0 - b) / (1.0 + b)
        c = kappa * x0 + (p - 1.0) * np.log(1.0 - x0**2)

        samples = np.zeros((N, 1))
        # Loop over number of samples
        for i in range(N):
            while True:
                Z = rng.beta((p - 1.0) / 2.0, (p - 1.0) / 2.0)
                U = rng.uniform(low=0.0, high=1.0)
                W = (1.0 - (1.0 + b) * Z) / (1.0 - (1.0 - b) * Z)

                if kappa * W + (p - 1.0) * np.log(1.0 - x0 * W) - c >= np.log(U):
                    samples[i] = W
                    break
        return samples

    def rand_von_mises_fisher(mu, kappa, N):
        p = len(mu)  # Dimension p
        mu = np.reshape(mu, (p, 1))  # Make sure that mu has a shape of px1
        samples = np.zeros((N, p))  # Array to store samples
        t = rand_t_marginal(kappa, p, N)  #  Component in the direction of mu (Nx1)
        xi = rand_uniform_hypersphere(N, p - 1)  # Component orthogonal to mu (Nx(p-1))

        # von-Mises-Fisher samples Nxp

        # Component in the direction of mu (Nx1). Note that here we are choosing an
        # intermediate mu = [1, 0, 0, 0, ..., 0] later we rotate to the desired mu below
        samples[:, [0]] = t
        # Component orthogonal to mu (Nx(p-1))
        samples[:, 1:] = (
            np.matlib.repmat(np.sqrt(1 - t**2), 1, p - 1) * xi
        )  # could probably be simplified to = np.sqrt(1 - t**2) * xi

        # Rotation of samples to desired mu

        O = null_space(mu.T)

        R = np.concatenate((mu, O), axis=1)
        samples = np.dot(R, samples.T).T
        return samples

    tiltvec = rand_von_mises_fisher(mu, kappa=kappa, N=2 * nNPCs)

    tiltnucv = tiltvec[:nNPCs]
    tiltcytv = tiltvec[nNPCs:]

    return tiltnucv, tiltcytv


def Multiple_coord(coords, symmet, nupIndex):  
    "Input is a list of coordinates seperated by NPC. output is all coordinates"

    nNPCs = len(coords)
    nRings = len(coords[0])

    coords_reshaped = np.reshape(coords, (nNPCs * nRings * symmet, 3))  # 3 Dim

    nupi = np.repeat(nupIndex, symmet)
    i = np.repeat(np.arange(nNPCs), symmet * nRings)


    return np.c_[(coords_reshaped, nupi, i)]


rotC = np.random.vonmises(0, 0, size=10000)  # 4*np.pi)
z = np.cos(rotC)  # z is x relative to on unit cicle. dist to NPC in z
x0 = np.sin(rotC)  # x is y. radius to NPC centre in x

# plt.scatter(x0, z, alpha = 0.01)
