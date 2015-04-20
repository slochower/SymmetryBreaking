#!/usr/bin/python

from math import *
from numpy import *
from scipy import ndimage
from timeit import default_timer as timer
import time

execfile('../common/figures.py')


def landscapes(x, shift=0):
    period = 2 * pi
    amplitude = 40
    window = 20

    # If the input is an array, then the energy function is evaluated at
    # every point in the array, and the complete landscape is returned.
    sawtooth = amplitude * ((x + shift) / period -
                            floor(1 / 2 + (x + shift) / period))
    padded = concatenate((sawtooth, sawtooth, sawtooth))
    smoothed = gaussian_filter(padded, window)
    if debug:
        plt.figure()
        plt.plot(range(len(padded)), padded, c='k', ls='--',
                 label='Original (padded)')
        plt.plot(range(len(smoothed)), smoothed, c='r', lw=2,
                 label='Smoothed')
        plt.legend()
        plt.title('Step function landscapes padded and smoothed')
        plt.show()
    return(smoothed[len(sawtooth):2*len(sawtooth)])


def gaussian_filter(x, N):
    return(ndimage.filters.gaussian_filter1d(input=x, sigma=N))


def find_nearest(array, value):
    idx = (abs(array - value)).argmin()
    # return array[idx]
    return idx


def force(i, dx, energy):
    force = -((energy[i+1] - energy[i]) / dx)
    if debug:
        print('Force = {}'.format(force))
    return force


def force_lookup(q, array, value):
    push = interp(value, q, array)

    if debug:
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        ax.plot(range(len(array[:, 1])), array[:, 1], color='g', lw=2,
                label='Force')
        ax.scatter(array[:, 0][nearest_index-1], array[:, 1][nearest_index-1],
                   color='k', s=40, zorder=10, alpha=0.5)
        ax.scatter(array[:, 0][nearest_index+1], array[:, 1][nearest_index+1],
                   color='k', s=40, zorder=10, alpha=0.5)
        ax.scatter(array[:, 0][nearest_index], push,
                   color='k', s=40, zorder=10)
        ax.set_title('Forces and interpolation')
        plt.show()

    return(push)

def energy_lookup(position, energy, value):
    light = interp(value, position, energy)
    return(light)


def simulate(x, dx, D, kT, dt, shift, forces, energy, i, **kwargs):

    ##################################
    # BROWNIAN MOTION ON THE LANDSCAPE
    ##################################
    positions = []        # List of walker positions
    jumps = []            # List of walker steps
    fluxes = []           # List of PBC transitions (+/- 1)
    # positions.append(x)   # Set initial position
    reason = -1           # No jump
    if debug:
        print('Starting walker at x = {}'.format(x))
    for t in range(timesteps):
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
        if (t % MC_time == 0 and t > MC_time and MC == True):
            E_transition = energy_lookup(q, energy[(i+1) % 2], x)
            E_now = energy_lookup(q, energy[i], x)
            # This might not be necessary:
            # E_step = energy_lookup(q, energy[i], candidate_step)
            if debug:
                print 'Current energy = {}'.format(E_now)
                print 'Energy step to new landscape = {}'.format(E_transition)
            delta = E_transition - E_now
            p_accept = exp(-delta/kT)
            r = random.random()
            if (delta < 0):
                print('Move because delta < 0 after {} steps.'.format(t))
                reason = 0
                break  # Kills the for() loop
            if (p_accept > r):
                print('Move because p_accept > r after {} steps'.format(t))
                reason = 1
                break  # Kills the for() loop
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(q[:-1], forces[i], x)
        new_x = x + (D / kT) * F * dt + g
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
            print 'Jumps are too big.'
        jumps.append(jump)
        x = new_x
    fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))
    ax = axes[0][0]
    ax.plot(q, energy[i], color='k', lw=2)
    ax.set_title('Energy landscape {} \n'
                 '{} steps taken, final $x = {:3.2f}$'.format(i, len(positions), positions[-1]), fontsize=fs['title'])
    if len(positions) > 0:
        colors = range(len(positions))
        if reason == 0:
            lbl = 'Jumped because $\Delta E < 0$'
        elif reason == 1:
            lbl = 'Jumped because $p_\\text{accept} < r$'
        else:
            lbl = 'No Monte Carlo jump'
        cax = ax.scatter(positions,
                         [energy_lookup(q, energy[i], positions[k])
                          for k in range(len(positions))],
                         c=colors, cmap=cm.jet, s=200, lw=0,
                         zorder=10, alpha=0.05,
                         label=lbl)
        # cbar = fig.colorbar(cax, ticks=[0, len(positions)])
        # cbar.ax.set_yticklabels(['Start', 'End'])
        # cbar.set_alpha(1)
        # cbar.draw_all()
        ax.legend(fontsize=fs['title'])
        ax.set_ylabel('Energy', fontsize=fs['axlabel'])
        ax.set_xlabel('Reaction coordinate', fontsize=fs['axlabel'])
        ax.tick_params(labelsize=fs['ticks'])
        ax.grid()
    plt.savefig(time.strftime("%Y-%m-%d")+'_landscape_{}.png'.format(i))
    plt.close()

    return(positions, jumps, fluxes)


#################################
# PARAMETERS
#################################
debug, deterministic, MC = False, False, True                # Tests
flashing = False
plot_together = True
interpolation_tests = False

num_landscapes = 20
# Number of protein states
# This is now a misnomer; it is really number of instances of 2 states.
# Should be even.
dx = 0.01
# Reaction coordinate spacing,
# although the walker is basically
# a continuous variable.
q = arange(0, 2 * pi, dx)
# Define reaction coordinate
dt = 1
# Time resolution (not really used)
timesteps = 1000000
# Number of steps for walker on each landscape
D = 0.01
# Arbitrary -- these two work together!
kT = 10
# Arbitrary -- these two work together!
MC_time = 10000
# Number of steps between Monte Carlo iterations

#################################
# SIMULATE
#################################

walker = [[] for i in range(num_landscapes)]
jumps = [[] for i in range(num_landscapes)]
net_flux = [[] for i in range(num_landscapes)]

shift = [0 if i % 2 == 0 else pi for i in range(num_landscapes)]
energy = [landscapes(q, shift[i]) for i in range(num_landscapes)]
if flashing:
    energy[1::2] = [array([mean(energy[0]) for i in range(len(q))]) for k in range(num_landscapes/2)]
forces = [[force(k, dx, energy[i]) for k in range(len(q)-1)]
          for i in range(num_landscapes)]
boltzmann = [exp(-energy[i]/kT) for i in range(num_landscapes)]
pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(num_landscapes)]

start = timer()
for i in range(num_landscapes):
    if i == 0:
        position = pi
    else:
        position = walker[i-1][-1]
    print('Walking on landscape {}'.format(i))
    tmp = simulate(position, dx, D, kT, dt, shift[i], forces, energy, i)

    print 'Adding {} to walker {}'.format(len(tmp[0]), i % 2)
    walker[i] = tmp[0]
    jumps[i] = tmp[1]
    net_flux[i] = tmp[2]

simulation_time = timer() - start
print('Simulation took {} seconds'.format(simulation_time))

fluxes = [sum(net_flux[i]) for i in range(len(net_flux))]
################################
# CONDENSE THE DATA
################################
print('Overall net flux = {}'.format(sum(fluxes)))
print('Ratio of steps on energy landscape A to B = {}'.
      format(float(len(hstack(walker[0::2]))) /
             float(len(hstack(walker[1::2])))))


#################################
# PLOT
#################################

if plot_together:
    fig, gs, axes = generate_axes_pad(nrows=2, ncols=2, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))
    for i in range(2):
        c = clrs[i % 9]
        ax = axes[i][0]
        ax.plot(q, energy[i], color=c, lw=2)
        ax.set_title('Energy landscape')

        ax = axes[i][1]
        # i = 0, then walker[0], walker[2], walker[4]
        # i = 1, then walker[1], walker[3], walker[5]
        counts, edges = histogram(hstack(walker[i::2]), range=(min(q), max(q)),
                                  bins=len(q), density=True)
        mids = (edges[1:] + edges[:-1]) / 2.
        ax.bar(mids, counts, color=c, edgecolor='none',
               width=mids[1] - mids[0], label='{} steps'.format(len(hstack(walker[i::2]))))
        ax.plot(q, pdf[i], color='k', lw=2)
        ax.set_title('Net flux = {}, kT = {}'.format(fluxes[i], kT))

        ax.legend()
    plt.show()


if interpolation_tests:

    from mpl_toolkits.mplot3d.axes3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X = arange(0, 2*pi, 0.25)
    Y = arange(0, 40, 1)
    x, y = meshgrid(X, Y)
    z = landscapes(x)
    z[-1] = landscapes(x[-1], shift=pi)
    z[1:-1] = [mean(energy[0]) for k in range(shape(z)[1])]
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    from scipy.interpolate import griddata
    grid_x, grid_y = mgrid[0:2*pi:0.25, 0:41:1]

    zero = array([q, [0]*629, landscapes(q)]).T
    one = array([q, [40]*629, landscapes(q, shift=pi)]).T
    values = vstack((zero, one))
    points = tuple([tuple(row) for row in values[:,0:2]])
    grid_z0 = griddata(array(points), values[:,2], (grid_x, grid_y), method='nearest')

    surf = ax.plot_surface(grid_x, grid_y, grid_z0, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    plt.show()    

