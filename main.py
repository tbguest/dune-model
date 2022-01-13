import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import random

import time

t0 = time.time()

"""
Assumptions: 
- stablitity of the prior bed state, allowing only the neighbours of the cell being 
acted upon to be checked
- wind blows to the north (+y)

Think carefully about atan(y/x), and which to apply the aspect ratio to
"""


def compute_slope(dh, dd, ar):
    """
    dh: y length [difference in no. of slabs]
    dd: x length [lattice cells]
    ar: aspect ratio
    """

    slope = np.arctan(dh / dd * ar) * 180 / np.pi

    return slope


def enforce_angle_of_repose(h, ix, iy, mode="add"):
    """
    ix, iy: idices of cell that has been added to or removed from
    Has a slab been removed or added?

    mode = "add" or "remove"
    """

    nn = {}
    nnSlopes = {}

    # this is elegant for readability, but makes the following logic more awkward...
    nn["N"] = [ix, iy + 1]
    nn["NE"] = [ix + 1, iy + 1]
    nn["E"] = [ix + 1, iy]
    nn["SE"] = [ix + 1, iy - 1]
    nn["S"] = [ix, iy - 1]
    nn["SW"] = [ix - 1, iy - 1]
    nn["W"] = [ix - 1, iy]
    nn["NW"] = [ix - 1, iy + 1]

    # apply wrap at boundaries
    for nbr in nn:
        if nn[nbr][0] < 0:
            nn[nbr][0] = nx - 1
        elif nn[nbr][0] > nx - 1:
            nn[nbr][0] = 0
        if nn[nbr][1] < 0:
            nn[nbr][1] = ny - 1
        elif nn[nbr][1] > ny - 1:
            nn[nbr][1] = 0

    # compute elevation difference for each neighbouring cell
    for nbr in nn:
        dh = h[ix, iy] - h[nn[nbr][0], nn[nbr][1]]
        nnSlopes[nbr] = compute_slope(dh, 1, slabAspectRatio)

    maxSlope = max(nnSlopes.values())

    stable = True
    if maxSlope > reposeAngle:

        # are there repeat maxima to choose from?
        candidates = []
        for nbr in nnSlopes:
            if nnSlopes[nbr] == maxSlope:
                candidates.append(nbr)

        # if more than one max, choose a fall direction at random
        iFall = int(np.floor(random.rand() * len(candidates)))
        # for the vanishingly small chance that rand() returns 1.0
        if iFall == len(candidates):
            iFall = iFall - 1
        ixFall = nn[candidates[iFall]][0]
        iyFall = nn[candidates[iFall]][1]

        if mode == "add":
            h[ix, iy] -= 1
            h[ixFall, iyFall] += 1
        elif mode == "remove":
            # NOTE that if this case is called, we've already ensured than h > 0
            h[ix, iy] += 1
            h[ixFall, iyFall] -= 1
        else:
            print("Invalid mode passed to enforce_angle_of_repose().")
            sys.exit()

        ix, iy = ixFall, iyFall
        stable = False

    return h, ix, iy, stable


# define the lattice
xmax = 1000
ymax = 1000
x = np.linspace(0, xmax - 1, xmax)
y = np.linspace(0, ymax - 1, ymax)
nx = len(x)
ny = len(y)

h = np.zeros((nx, ny))

# "slab" aspect ratio (< of repose = atan(2/3))
slabAspectRatio = 1.0 / 3.0
meanSlabs = 3  # mean number of slabs/lattice point

Ps = 0.6
Pns = 0.4
L = 5

reposeAngle = 30  # degrees
shadowAngle = 15  # degrees

tmax = 500

# fig, ax = plt.subplots(nrows=1, ncols=1)


# Setup step:

# populate lattice with randomly placed slabs
for _ in range(nx * ny * meanSlabs):

    # pick a lattice point at random
    ix = int(np.round((nx - 1) * random.rand()))
    iy = int(np.round((ny - 1) * random.rand()))

    # add a "slab"
    h[ix, iy] = h[ix, iy] + 1

    # check the angle of repose
    stable = False
    while not stable:
        h, ix, iy, stable = enforce_angle_of_repose(h, ix, iy, mode="add")


# plt.plot(x, h)
# plt.contour(h)
# plt.show()

######## main model ########

tcount = 0
allcount = 0


# def check_shadow_zone():

#     return


while tcount < tmax:

    # saltation step:

    # pick an index at random
    ix = int(np.round((nx - 1) * random.rand()))
    iy = int(np.round((ny - 1) * random.rand()))

    # if the substrate is not exposed:
    if h[ix, iy] > 0:

        # check if in shadow zone:

        """
        - check all upwind coords (all coords on lattice) to see if threshold angle is met
        - associate an index with every upwind cell
        - compute distances and relative heights for each index

        should omit self to avoid / by 0
        dUpwind[ix]

        """

        # compute the upwind distance and angles to each lattice cell
        dUpwind = (iy - y) % nx
        thetas = np.arctan((h - h[ix, iy]) / dUpwind * slabAspectRatio) * 180 / np.pi

        if np.nansum(thetas > shadowAngle) == 0:
            # not in shadow zone; proceed

            # select grain to saltate
            h[ix, iy] = h[ix, iy] - 1

            # retain original indices, since they may be mutated in repose step
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
                if iy > ny - 1:
                    iy = iy % ny

                # compute the upwind distance and angles to each lattice cell
                dUpwind = (iy - y) % ny
                thetas = (
                    np.arctan((h - h[ix, iy]) / dUpwind * slabAspectRatio) * 180 / np.pi
                )

                if np.nansum(thetas > shadowAngle) > 0:
                    h[ix, iy] += 1  # deposited with P=1
                    transport = False

                # if there's sand at the condidate deposition site
                elif h[ix, iy] > 0:

                    # slab deposited with P=Ps
                    if random.rand() < Ps:
                        h[ix, iy] = h[ix, iy] + 1
                        transport = False
                else:
                    # slab deposited with P=Pns
                    if random.rand() < Pns:
                        h[ix, iy] = h[ix, iy] + 1
                        transport = False

            # check the angle of repose
            stable = False
            while not stable:
                h, ix, iy, stable = enforce_angle_of_repose(h, ix, iy, mode="add")

    tcount = tcount + 1 / nx
    allcount = allcount + 1

    # if allcount % 5000 == 1:
    #     ax.plot(x, h[] + tcount / 20 - 1, "k", "linewidth", 1.5)


print((time.time() - t0) / 60, "mins")

# plt.plot(y, h[0, :])
plt.contour(h)
plt.show()


# # probability of deposition

# # check if the slab landed in a shadow zone
# Icheck = ix - np.arange(1, nx - 2, step=1)
# iDownwind = np.argwhere(Icheck < 1)
# for trsh in range(len(iDownwind)):
#     Icheck[iDownwind[trsh]] = Icheck[iDownwind[trsh]] + nx - 1

# counter = 0
# for jj in range(len(Icheck)):
#     if (
#         np.arcatan2(jj * 3 / 2, h[Icheck[jj]]) > 15 * 180 / np.pi
#     ):  # (h(Icheck(jj))/(jj*3/2)) > 15*180/pi
#         counter = counter + 1

# if counter > 0:
#     h[ix] = h[ix] + 1  # deposited with P=1
#     transport = False

# elif h[ix] > 0:
#     if random.rand() < Ps:
#         h[ix] = h[ix] + 1
#         transport = False
# else:
#     if random.rand() < Pns:
#         h[ix] = h[ix] + 1
#         transport = False

