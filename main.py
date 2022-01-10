import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import random


"""
Assume: stablitity of the prior bed state, allowing only the neighbours of the cell being 
acted upon to be checked

Think carefully about atan(y/x), and which to apply the aspect ratio to
"""


def enforce_angle_of_repose(h, x, ix, iy):

    """Has a slab been removed or added?"""

    # handle boundary cases for nearest neighbour (nn) cells
    if ix == 0 and iy == 0:
        """Not sure this is the most effective way to do this"""
        nnN = 1
        nnNE = 
        nnE = 1
        nnSE = 
        nnS = ny - 1
        nnSW = 
        nnW = nx - 1
        nnNW = 
    elif ix == nx - 1:
        nnW = nx - 2
        nnE = 0

    if ix == 0:
        nnW = nx - 1
        nnE = 1
    elif ix == nx - 1:
        nnW = nx - 2
        nnE = 0

    if iy == 0:
        nnS = ny - 1
        nnN = 1
    elif iy == ny - 1:
        nnS = ny - 2
        nnN = 0

    else:
        Ilo = ix - 1
        Ihi = ix + 1

    return h

# define the lattice
xmax = 1000
ymax= 1000
x = np.linspace(0, xmax - 1, xmax)
y = np.linspace(0, ymax - 1, ymax)
nx = len(x)
ny = len(y)

h = np.zeros(nx,)

# "slab" aspect ratio (< of repose = atan(2/3))
slabAspectRatio = 1.0 / 3.0
meanSlabs = 1  # mean number of slabs/lattice point

Ps = 0.6
Pns = 0.4
L = 5.0

reposeAngle = 30  # degrees
shadowAngle = 15  # degrees

tmax = 10

fig, ax = plt.subplots(nrows=1, ncols=1)


# Setup step:

# populate lattice with randomly placed slabs
for ii, _ in enumerate(np.repeat(x, meanSlabs)):

    # pick a lattice point at random
    ix = int(np.round((nx - 1) * random.rand()))

    # add a "slab"
    h[ix] = h[ix] + 1

    # check the angle of repose
    stable = False  # conditions not yet met
    while not stable:
        if ix == 0:
            Ilo = nx - 1
            Ihi = 1
        elif ix == nx - 1:
            Ilo = nx - 2
            Ihi = 0
        else:
            Ilo = ix - 1
            Ihi = ix + 1

        d1 = h[ix] - h[Ilo]
        d2 = h[ix] - h[Ihi]

        if np.abs(d1) > 1 and np.abs(d2) > 1:
            h[ix] = h[ix] - 1
            if d1 == d2:
                Ipickone = np.round(random.rand())
                if Ipickone == 1.0:
                    h[Ilo] = h[Ilo] + 1
                    ix = Ilo
                else:
                    h[Ihi] = h[Ihi] + 1
                    ix = Ihi

        elif np.abs(d1) > 1 or np.abs(d2) > 1:
            h[ix] = h[ix] - 1
            if np.abs(d1) > np.abs(d2):
                h[Ilo] = h[Ilo] + 1
                ix = Ilo
            else:
                h[Ihi] = h[Ihi] + 1
                ix = Ihi

        # condition has been met
        else:
            stable = True


plt.plot(x, h)

plt.show()


######## main model ########

tcount = 0
allcount = 0

while tcount < tmax:

    print(tcount)

    # saltation step:

    # pick an index at random
    ix = int(np.round((nx - 1) * random.rand()))

    # if the substrate is not exposed:
    if h[ix] > 0:

        # check if in shadow zone:

        """
        - check all upwind coords (all coords on lattice) to see if threshold angle is met
        - associate an index with every upwind cell
        - compute distances and relative heights for each index
        
        should omit self to avoid / by 0
        dUpwind[ix]
         
        """

        # compute the upwind distance and angles to each lattice cell
        dUpwind = (ix - x) % nx
        thetas = np.arctan((h - h[ix]) / dUpwind * slabAspectRatio) * 180 / np.pi

        if np.sum(thetas > shadowAngle) == 0:
            # not in shadow zone; proceed

            # select grain to saltate
            h[ix] = h[ix] - 1

            # check angle of repose:
            stable = False  # conditions not yet met
            while not stable:

                # wrap conditions
                if ix == 0:
                    Ilo = nx - 1
                    Ihi = 1
                elif ix == nx - 1:
                    Ilo = nx - 2
                    Ihi = 0
                else:
                    Ilo = ix - 1
                    Ihi = ix + 1

                d1 = h[ix] - h[Ilo]
                d2 = h[ix] - h[Ihi]
                if np.abs(d1) > 1 and np.abs(d2) > 1:
                    h[ix] = h[ix] + 1
                    if d1 == d2:
                        Ipickone = np.round(1 + (np.random.rand()))

                        # TODO: I don't see that interative instabilities that get introduces get properly handled here..

                        if Ipickone == 1.0:
                            h[Ilo] = h[Ilo] - 1
                            ix = Ilo
                        else:
                            h[Ihi] = h[Ihi] - 1
                            ix = Ihi

                elif np.abs(d1) > 1 or np.abs(d2) > 1:
                    h[ix] = h[ix] + 1
                    if np.abs(d1) > np.abs(d2):
                        h[Ilo] = h[Ilo] - 1
                        ix = Ilo
                    else:
                        h[Ihi] = h[Ihi] - 1
                        ix = Ihi

                # condition has been met
                else:
                    stable = True

            # grain saltates:
            transport = True

            """Wait! ix may have been mutated above!"""
            while transport:
                ix = ix + L
                if ix > nx - 1:
                    ix = ix % nx

                # prob of deposition

                # fixst, check if in shadow zone
                Icheck = ix - np.arange(1, nx - 2, step=1)
                iDownwind = np.argwhere(Icheck < 1)
                for trsh in range(len(iDownwind)):
                    Icheck[iDownwind[trsh]] = Icheck[iDownwind[trsh]] + nx - 1

                counter = 0
                for jj in range(len(Icheck)):
                    if (
                        np.arcatan2(jj * 3 / 2, h[Icheck[jj]]) > 15 * 180 / np.pi
                    ):  # (h(Icheck(jj))/(jj*3/2)) > 15*180/pi
                        counter = counter + 1

                if counter > 0:
                    h[ix] = h[ix] + 1  # deposited with P=1
                    transport = False

                elif h[ix] > 0:
                    if random.rand() < Ps:
                        h[ix] = h[ix] + 1
                        transport = False
                else:
                    if random.rand() < Pns:
                        h[ix] = h[ix] + 1
                        transport = False

            # check angle of repose
            stable = False  # conditions not yet met
            while not stable:
                if ix == 0:
                    Ilo = nx - 1
                    Ihi = 1
                elif ix == nx - 1:
                    Ilo = nx - 2
                    Ihi = 0
                else:
                    Ilo = ix - 1
                    Ihi = ix + 1

                d1 = h[ix] - h[Ilo]
                d2 = h[ix] - h[Ihi]
                if np.abs(d1) > 1 and np.abs(d2) > 1:
                    h[ix] = h[ix] - 1
                    if d1 == d2:
                        Ipickone = np.round(1 + random.rand())
                        if Ipickone == 1.0:
                            h[Ilo] = h[Ilo] + 1
                            ix = Ilo
                        else:
                            h[Ihi] = h[Ihi] + 1
                            ix = Ihi

                elif np.abs(d1) > 1 or np.abs(d2) > 1:
                    h[ix] = h[ix] - 1
                    if np.abs(d1) > np.abs(d2):
                        h[Ilo] = h[Ilo] + 1
                        ix = Ilo
                    else:
                        h[Ihi] = h[Ihi] + 1
                        ix = Ihi

                # condition has been met
                else:
                    stable = True

    tcount = tcount + 1 / nx
    allcount = allcount + 1

    if np.abs(np.max(np.diff(h.flatten()))) > 1:
        print("repose condition not met")
        break

    ind = 0
    if allcount % 5000 == 1:
        ax.plot(x, h + tcount / 20 - 1, "k", "linewidth", 1.5)


plt.plot(x, h)

plt.show()
