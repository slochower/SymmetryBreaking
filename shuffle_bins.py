#!/usr/bin/python

from math import *
from numpy import *

from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory
from timeit import default_timer as timer

execfile('../common/figures.py')
'''
def landscapes(x, shift=0):
    period = 1.8 * pi
    amplitude = 10
    window = 10

    # If the input is an array, then the energy function is evaluated at
    # every point in the array, and the complete landscape is returned.
    if isinstance(x, ndarray):
        sawtooth = amplitude * ((x + shift) / period -
                                floor(1 / 2 + (x + shift) / period))
        # This is necessary because the running average will drop the 
        # first window-1 points, so we take advantage of PBCs.
        tmp = concatenate((sawtooth[-window:], sawtooth))
        smoothed = pd.rolling_mean(tmp, window=window)
        if debug:
            print('landscapes() called with array.')
        return(smoothed[window:])
    # If the input is a point (e.g., the walker is at some position and
    # and we need to spit out the energy at that specific position) then
    # we have to generate the landscape a little before and a little after
    # that position, so it can be smoothed.
    if isinstance(x, float):
        # xx = arange(x - window * dx, x + window * dx, dx)
        # The results are not consistent with arange().
        xx = linspace(x - window * dx, x + window * dx, 2*window+1)
        sawtooth = amplitude * ((xx + shift) / period -
                                floor(1 / 2 + (xx + shift) / period))
        smoothed = pd.rolling_mean(sawtooth, window=10)
        if debug:
            print('landscapes() called with float.')
            print('x = {}'.format(x))

        return(sawtooth[window:-window])
    else:
        print('Erorr computing energy landscape.')
        print('type(x) = {}'.format(type(x)))
'''

def landscapes(x, shift=0):
    period = 2 * pi
    amplitude = 10
    if debug:
        print('landscapes() reporting in with x = {}'.format(x))
    sawtooth = amplitude * ((x + shift) / period -
                            floor(1 / 2 + (x + shift) / period))
    if debug:
        print('landscapes() reporting in with E = {}'.format(sawtooth))
    return(sawtooth)


def landscapes(x, shift=0, noise=0):
    return ((cos(2 * x + shift * pi))/2)


def find_nearest(array, value):
    idx = (abs(array - value)).argmin()
    return array[idx]


def force(x, dx, s):
    if debug:
        print('Calculating force with x, x+dx, shift = {}, {}, {}'.
              format(x, x+dx, s))
    force = -(landscapes(x + dx, shift=s) -
              landscapes(x, shift=s)) / dx
    if debug:
        print('Force = {}'.format(force))
    return float(force)


def smooth(x, N):
    # Boundary effects make this not so ideal.
    return convolve(x, ones((N,)) / N, mode='valid')


def smooth_average(x, N):
    y = zeros((len(x),))
    for ctr in range(len(x)):
        y[ctr] = sum(x[ctr:(ctr + N)])
    return y / N


def simulate(x, dx, D, kT, dt, shift, **kwargs):

    ##################################
    # BROWNIAN MOTION ON THE LANDSCAPE
    ##################################
    positions = []        # List of walker positions
    jumps = []            # List of walker steps
    fluxes = []           # List of PBC transitions (+/- 1)
    forces = []           # List of forces on walker
    # positions.append(x)   # Set initial position

    if debug:
        print('Starting walker at x = {}'.format(x))
    for t in range(timesteps):
        swap = 0
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force(x, dx, shift)
        new_x = x + (D / kT) * F * dt + g
        if deterministic:
            new_x = x + F * dt
        if debug:
            print 't = {} \t x = {}, force * dt = {}, new_x = {}' \
                  .format(float(t), float(x),
                          float(F * dt), float(new_x))

        if new_x > max(q):
            new_x = min(q) + (new_x - max(q))
            swap = 1
            if debug:
                print 'Wrapping {} to {}' \
                      .format(new_x, min(q) + (new_x - max(q)))
                print 'Setting swap to {}'.format(swap)
                print 'Jump = {}'.format(new_x - x + swap * (max(q) - min(q)))
        elif new_x < min(q):
            if debug:
                print 'Wrapping {} to {}'.format(new_x, max(q) -
                                                 abs((new_x - min(q))))
            new_x = max(q) - abs(new_x - min(q))
            swap = -1
        else:
            swap = 0
        fluxes.append(swap)
        positions.append(new_x)
        jump = new_x - x + swap * (max(q) - min(q))
        if (abs(jump) > max(q)):
            print new_x
            print x
            print swap
            print 'Jump =  {} '.format(new_x - x + swap * (max(q) - min(q)))
            print 'WARNING: JUMPED ACROSS AN ENTIRE LANDSCAPE.'
        jumps.append(jump)
        forces.append(F)
        x = new_x

    if debug:
        print 'Final x = {}'.format(x)
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        counts, edges = histogram(jumps, bins=100, density=True)
        mids = (edges[1:] + edges[:-1]) / 2.
        ax.bar(mids, counts, color='g', edgecolor='none',
               width=mids[1] - mids[0])
        ax.set_title('Steps')
        plt.savefig('figure.png')
        plt.close()

    return(positions, jumps, fluxes, forces)


#################################
# PARAMETERS
#################################
debug, deterministic = False, False                 # Tests
plot_together = True

num_landscapes = 3
# Number of protein states
dx = 0.01
# Reaction coordinate spacing,
# although the walker is basically
# a continuous variable.
q = arange(0, 2 * pi, dx)
# Define reaction coordinate
dt = 1
# Time resolution (not really used)
timesteps = 100000
# Number of steps for walker on each landscape
D = 0.01
# Arbitrary -- these two work together!
kT = 10
# Arbitrary -- these two work together!


walker = zeros((num_landscapes, timesteps))
jumps = zeros((num_landscapes, timesteps))
net_flux = zeros((num_landscapes, timesteps))
forces = zeros((num_landscapes, timesteps))

shift = [0 if i % 2 == 0 else pi for i in range(num_landscapes)]
energy = [landscapes(q, shift[i]) for i in range(num_landscapes)]
boltzmann = [exp(-energy[i]/kT) for i in range(num_landscapes)]
pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(num_landscapes)]

start = timer()
for i in range(num_landscapes):
    if i == 0:
        position = pi
    else:
        position = walker[i-1][-1]
    print('Walking on landscape {}'.format(i))
    walker[i], jumps[i], net_flux[i], forces[i] = simulate(position, dx,
                                                 D, kT, dt, shift[i])
simulation_time = timer() - start
print('Simulation took {} seconds'.format(simulation_time))

fluxes = sum(net_flux, axis=1)
print('Overall net flux = {}'.format(sum(net_flux)))

if plot_together:
    fig, gs, axes = generate_axes_pad(nrows=num_landscapes, ncols=3, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))

    ##################################
    # PLOTTING LANDSCAPES TOGETHER
    ##################################
    for i in range(num_landscapes):
        c = clrs[i % 9]
        ax = axes[i][0]
        ax.scatter(walker[i][0], landscapes(walker[i][0], shift=shift[i]),
                   color='k', s=40, zorder=10, label='Start')
        ax.scatter(walker[i][-1], landscapes(walker[i][-1], shift=shift[i]),
                   color='k', alpha=0.5, s=40, zorder=10, label='End')
        ax.plot(q, energy[i], color=c, lw=2)
        ax.set_title('Energy landscape')

        ax = axes[i][1]
        counts, edges = histogram(walker[i], range=(min(q), max(q)),
                                  bins=len(q), density=True)
        mids = (edges[1:] + edges[:-1]) / 2.
        ax.bar(mids, counts, color=c, edgecolor='none',
               width=mids[1] - mids[0])
        ax.plot(q, pdf[i], color='k', lw=2, label='kT = {}'.format(kT))
        ax.set_title('Net flux = {}, kT = {}'.format(fluxes[i], kT))

        ax = axes[i][2]
        counts, edges = histogram(jumps[i], bins=100, density=True)
        mids = (edges[1:] + edges[:-1]) / 2.
        ax.bar(mids, counts, color=c, edgecolor='none',
               width=mids[1] - mids[0])
        ax.set_title('Step sizes')
        ax.legend()
plt.show()
