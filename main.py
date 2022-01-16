"""
Werner's (1995) dune model, with additions from Momiji et al. (2000).

Dunes are built from piled "slabs" representing sand on a one- or two-dimensional lattice, 
whose edges are connected with periodic boundary conditions. The number of sand slabs is proportional
to the surface height, and their movement comprises one-directional transport together with avalanching dynamics.
A sand slab is individually and randomly chosen for transport (erosion) from among all the slabs on the
surface. The slab is moved a specific number of lattice sites, thus introducing a transport length (L), and is
deposited at a site with a probability P_s. If the slab is not deposited, then it is
repeatedly moved over L sites until deposition. Subsequently another slab is chosen randomly for transport.
Shadow zones are introduced in the lee of dunes. If a slab is transported to a site in a shadow zone, it is
deposited with unit probability. This process is repeated to construct the time evolution of the dune-field
surface.

Assumptions: 
- Stablitity of the prior bed state, allowing only the neighbours of the cell being acted upon to be checked
- Wind blows to the north (positive y direction)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import time

t0 = time.time()


### KNOBS ###
LATTICE_X = 1  # lattice dimensions
LATTICE_Y = 1000  # lattice dimensions
SAR = 1 / 3  # slab aspect ratio (y/x)
MEAN_SLAB_HEIGHT = 10  # mean number of slabs/lattice point
P_SAND = 0.6  # probability of slab deposition atop another sand slab
P_NOSAND = 0.4  # probability of slab deposition on bare substrate
L = 5  # slab saltation length
REPOSE_ANGLE = np.arctan(2 / 3) * 180 / np.pi  # angle of repose for sand (ca. 33 deg)
SHADOW_ANGLE = 15  # [deg] angle defining the shadow zone
TIMESTEPS = 150


def compute_slopes(dh, dd, ar):
    """
    Args:
        dh (1xn array): height difference [number of slabs]
        dd (1xn array): horizonral separation [number of lattice cells]
        ar (float): aspect ratio (y/x)
    """
    slopes = np.arctan2(dh * ar, dd) * 180 / np.pi
    return slopes


def enforce_angle_of_repose(h, ix, iy, mode="add"):
    """Ensure slope stability threshold is not exceeded after adding or removing a slab.

    Args:
        h (array): height of sand slabs at all lattice cells.
        ix, iy (int): indices of cell that has been added to or removed from.
        mode (str): "add" or "remove"

    Returns:
        h, ix, iy: updated values if threshold was exceeded.
        stable (boolean): True if no slab movements were required. 
    """

    # define indices of neareast neighbour (nn) cells, defined clockwise beginning with compass north (+y)
    nn = np.zeros((2, 8)).astype(int)
    nn[:, 0] = ix, iy + 1
    nn[:, 1] = ix + 1, iy + 1
    nn[:, 2] = ix + 1, iy
    nn[:, 3] = ix + 1, iy - 1
    nn[:, 4] = ix, iy - 1
    nn[:, 5] = ix - 1, iy - 1
    nn[:, 6] = ix - 1, iy
    nn[:, 7] = ix - 1, iy + 1

    # apply wrap at boundaries
    for ind, _ in enumerate(nn[0, :]):
        if nn[0, ind] < 0:
            nn[0, ind] = LATTICE_X - 1
        elif nn[0, ind] > LATTICE_X - 1:
            nn[0, ind] = 0
        if nn[1, ind] < 0:
            nn[1, ind] = LATTICE_Y - 1
        elif nn[1, ind] > LATTICE_Y - 1:
            nn[1, ind] = 0

    # compute elevation differences and slopes for each neighbouring cell
    dh = h[ix, iy] - h[nn[0, :], nn[1, :]]
    nnSlopes = compute_slopes(dh, np.ones(nn.shape[1]), SAR)
    maxSlopeIndex = np.argmax(np.abs(nnSlopes))
    maxSlope = nnSlopes[maxSlopeIndex]

    stable = True
    if np.abs(maxSlope) > REPOSE_ANGLE:

        # are there repeat maxima to choose from?
        isMax = maxSlope - nnSlopes == 0
        nMax = np.sum(isMax)
        if nMax > 1:

            # if more than one max, choose a fall direction at random
            iFall = int(np.floor(np.random.rand() * nMax))

            # wrap the index if rand() returns 1.0
            if iFall == nMax:
                iFall = 0
            isMaxIndex = np.where(isMax == True)[0]
            ixFall = nn[0, isMaxIndex[iFall]]
            iyFall = nn[1, isMaxIndex[iFall]]

        else:
            isMaxIndex = np.where(isMax == True)[0]
            ixFall = nn[0, isMaxIndex[0]]
            iyFall = nn[1, isMaxIndex[0]]

        if mode == "add":
            # check to ensure angle is positive
            if nnSlopes[isMaxIndex[0]] < 0:
                print("Encountered a negative angle on a slab removal step. Invalid.")
                sys.exit()

            h[ix, iy] -= 1
            h[ixFall, iyFall] += 1

        elif mode == "remove":
            # check to ensure angle is negative
            if nnSlopes[isMaxIndex[0]] > 0:
                print("Encountered a positive angle on a slab addition step. Invalid.")
                sys.exit()

            h[ix, iy] += 1
            h[ixFall, iyFall] -= 1

        else:
            print("Invalid mode passed to enforce_angle_of_repose().")
            sys.exit()

        ix, iy = ixFall, iyFall

        stable = False

    return h, ix, iy, stable


def initialize_model():
    """Place slabs at random until the chosen mean slab height is achieved.
    
    TODO: Optimize for larger mean heights"""

    h = np.zeros((LATTICE_X, LATTICE_Y))

    # populate lattice with randomly placed slabs
    for _ in range(LATTICE_X * LATTICE_Y * MEAN_SLAB_HEIGHT):

        # pick a lattice point at random
        ix = int(np.round((LATTICE_X - 1) * np.random.rand()))
        iy = int(np.round((LATTICE_Y - 1) * np.random.rand()))

        # add a "slab"
        h[ix, iy] = h[ix, iy] + 1

        # check the angle of repose
        stable = False
        while not stable:
            h, ix, iy, stable = enforce_angle_of_repose(h, ix, iy, mode="add")

    # h = MEAN_SLAB_HEIGHT * np.ones((LATTICE_Y, LATTICE_X))

    return h


if __name__ == "__main__":

    x = np.linspace(0, LATTICE_X - 1, LATTICE_X)
    y = np.linspace(0, LATTICE_Y - 1, LATTICE_Y)
    h = initialize_model()

    timestep = 0

    # for plotting time evolution in 1-d
    if LATTICE_X == 1 or LATTICE_Y == 1:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 3))

    while timestep < TIMESTEPS:

        # pick an index at random
        ix = int(np.round((LATTICE_X - 1) * np.random.rand()))
        iy = int(np.round((LATTICE_Y - 1) * np.random.rand()))

        # if the substrate is not exposed:
        if h[ix, iy] > 0:

            # compute the upwind distance and angles to each lattice cell
            dUpwind = (iy - y) % LATTICE_Y
            thetas = compute_slopes(h - h[ix, iy], dUpwind, SAR)

            # if not in the shadow zone
            if np.nansum(thetas > SHADOW_ANGLE) == 0:

                # remove the slab
                h[ix, iy] = h[ix, iy] - 1

                # retain original indices, since they may be mutated in the repose step
                ix0, iy0 = ix, iy

                # check the angle of repose
                stable = False
                while not stable:
                    h, ix0, iy0, stable = enforce_angle_of_repose(
                        h, ix0, iy0, mode="remove"
                    )

                # slab saltation step
                transport = True
                while transport:
                    iy = iy + L
                    if iy > LATTICE_Y - 1:
                        iy = iy % LATTICE_Y

                    # compute the upwind distance and angles to each lattice cell
                    dUpwind = (iy - y) % LATTICE_Y
                    thetas = compute_slopes(h - h[ix, iy], dUpwind, SAR)

                    # if in a shadow zone, deposit with unit probability (P=1)
                    if np.nansum(thetas > SHADOW_ANGLE) > 0:
                        h[ix, iy] += 1
                        transport = False

                    # if the cell contains one or more sand slabs, deposit with probability P=P_SAND
                    elif h[ix, iy] > 0:
                        if np.random.rand() < P_SAND:
                            h[ix, iy] = h[ix, iy] + 1
                            transport = False

                    # if bare substrate, deposit with probability P=P_NOSAND
                    else:
                        if np.random.rand() < P_NOSAND:
                            h[ix, iy] = h[ix, iy] + 1
                            transport = False

                # check the angle of repose
                stable = False
                while not stable:
                    h, ix, iy, stable = enforce_angle_of_repose(h, ix, iy, mode="add")

        timestep = timestep + 1 / (LATTICE_X * LATTICE_Y)

        # for plotting time evolution in 1-d
        if LATTICE_X == 1 or LATTICE_Y == 1:
            if timestep % 10 < 1 / (LATTICE_X * LATTICE_Y) or timestep == 1 / (
                LATTICE_X * LATTICE_Y
            ):
                ax.plot(y, h[0, :] + timestep - 1, "k", linewidth=2.0)

    ax.set_ylabel("y")
    ax.set_xlabel("x")
    ax.set_xlim([np.min(y), np.max(y)])
    fig.tight_layout()
    plt.savefig("dunes1d.png")

    print((time.time() - t0) / 60, "mins")

    fig2, ax2 = plt.subplots(nrows=1, ncols=1)
    if LATTICE_X == 1 or LATTICE_Y == 1:
        ax2.plot(y, h[0, :])
    else:
        ax2.contour(h)
    plt.show()
