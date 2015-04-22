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
    positions = []
    # List of walker positions on this landscape
    jumps = []
    # List of walker steps (could be derived from positions)
    fluxes = []
    # Crossings of the flux barrier
    # -1: left crossing of the boundary
    #  0: no crossing this step
    # +1: right crossing of the boundary
    reason = -1
    # Keep track of why MC attempted accepted
    # -1: step rejected
    #  0: delta E < 0
    # +1: p_accept > r

    if (i % 2) == 0: 
        state = 0
    else:
        state = 1


    if debug:
        print('Starting walker at x = {}'.format(x))
    for t in range(timesteps):
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
        if (MC is True and (t % MC_interval) == 0 and t > MC_interval):
            E_transition = energy_lookup(q, energy[(state+1) % 2], x)
            E_now = energy_lookup(q, energy[state], x)
            if debug:
                print 'Current energy = {}'.format(E_now)
                print 'Energy step to new landscape = {}'.format(E_transition)
            delta = E_transition - E_now
            p_accept = exp(-delta/kT)
            r = random.random()
            if (delta < 0):
                if debug:
                    print('Move because delta < 0 after {} steps.'.format(t))
                reason = 0
                break  # Kills the for() loop
            if (p_accept > r):
                if debug:
                    print('Move because p_accept > r after {} steps'.format(t))
                reason = 1
                break  # Kills the for() loop
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(q[:-1], forces[state], x)
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
        ####################
        # CLEANUP AND RECORD
        ####################
        fluxes.append(swap)
        positions.append(new_x)
        jump = new_x - x + swap * (max(q) - min(q))
        if (abs(jump) > max(q)):
            print 'Jumps are too big.'
        jumps.append(jump)
        x = new_x
    if plot_every is True:
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        ax.plot(q, energy[state], color='k', lw=2)
        ax.set_title('Energy landscape {} \n'
                     '{} steps taken, final $x = {:3.2f}$'.
                     format(state, len(positions), positions[-1]),
                     fontsize=fs['title'])
        if len(positions) > 0:
            colors = range(len(positions))
            if reason == 0:
                lbl = 'Jumped because $\Delta E < 0$'
            elif reason == 1:
                lbl = 'Jumped because $p_\\text{accept} < r$'
            else:
                lbl = 'No Monte Carlo jump'
            cax = ax.scatter(positions,
                             [energy_lookup(q, energy[state], positions[k])
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

    if (t == (timesteps - 1) and MC is True):
        print 'There were no transitions between landscapes during this'
        +' simulation. We are deviating from what we know and '
        +' understand.'
    if (t == (timesteps - 1) and MC is False):
        print 'Finished running simulation on landscape.'
        print 'This is the expected behavior.'
    # Should probably clean up the variables...
    return(positions, jumps, fluxes, t)


#################################
# PARAMETERS
#################################
debug, deterministic, MC = False, False, True
flashing = False
plot_together, plot_every = True, False
interpolation_tests = False

# num_landscapes = 20
# Number of protein states
# This is now a misnomer; it is really number of instances of 2 states.
# Should be even.
dx = 0.01
# Reaction coordinate spacing, although the walker is basically 
# a continuous variable.
q = arange(0, 2 * pi, dx)
# Define reaction coordinate
dt = 1
# Time resolution
# Don't change this because it is tied in with the array spacing
# e.g. positions[1] - positions[0] is the step from time 0 to time 1.
timesteps = 10000
# Maximum number of steps for walker on each landscape
total_timesteps = 1000000
# Cumulative number of steps for the simulation (across landscapes)
D = 0.01
# Arbitrary -- these two work together!
kT = 10
# Arbitrary -- these two work together!
MC_interval = 100
# Number of steps between Monte Carlo iterations

#################################
# SIMULATE
#################################

shift = [0, pi]
energy = [landscapes(q, shift[i]) for i in range(len(shift))]
if flashing:
    energy[1] = array([mean(energy[0])]*len(energy[0]))

forces = [[force(k, dx, energy[i]) for k in range(len(q)-1)]
          for i in range(len(energy))]
boltzmann = [exp(-energy[i]/kT) for i in range(len(energy))]
pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(len(boltzmann))]

start = timer()
#####################################################
# DEPRECATING LOOPING OVER A SET NUMBER OF LANDSCAPES
#####################################################
# for i in range(num_landscapes):
#     if i == 0:
#         position = pi
#     else:
#         position = walker[i-1][-1]
#     print('Walking on landscape {}'.format(i))
#     tmp = simulate(position, dx, D, kT, dt, shift[i], forces, energy, i)

#     print 'Adding {} to walker {}'.format(len(tmp[0]), i % 2)
#     walker[i] = tmp[0]
#     jumps[i] = tmp[1]
#     net_flux[i] = tmp[2]

# Maximum length on a single landscape, by maximum number of landscapes

walker = empty((total_timesteps/MC_interval, timesteps))
jumps = empty((total_timesteps/MC_interval, timesteps))
net_flux = empty((total_timesteps/MC_interval, timesteps))

i = 0
steps_executed = []
while True:
    if i == 0:
        position = pi
    elif i == 1:
        position = walker[-1][-1]
    else:
        position = walker[-1][-1]
        # This should pull out the last position of the last walker...
    print('Walking on landscape {}'.format(i))
    this_run = simulate(position, dx, D, kT, dt, shift, forces, energy, i)
    if i == 0:
        walker[0,0:len(this_run[0])] = this_run[0]
        jumps[0,0:len(this_run[1])] = this_run[1]
        net_flux[0,0:len(this_run[2])] = this_run[2]
        steps_executed.append(this_run[3])
    else:
    #    walker = vstack((walker, this_run[0]))
    #    jumps = vstack((jumps, this_run[1]))
    #    net_flux = vstack((net_flux, this_run[2]))
    #    import pdb; pdb.set_trace()
    #    walker.append(list(this_run[0]))
    #    jumps.append(list(this_run[1]))
    #    net_flux.append(list(this_run[2]))
    #    steps_executed.append(this_run[3])
        walker[i,0:len(this_run[0])] = this_run[0]
        jumps[i,0:len(this_run[1])] = this_run[1]
        net_flux[i,0:len(this_run[2])] = this_run[2]
        steps_executed.append(this_run[3])

    print('Total steps so far: {}'.format(sum(steps_executed)))

    if ( sum(steps_executed) >= total_timesteps):
        break
    i += 1


simulation_time = timer() - start
print('Simulation took {} seconds'.format(simulation_time))
print('Total landscapes sampled: {}'.format(i))


################################
# CONDENSE THE DATA
################################
# Eliminate all the zero rows... do this in new arrays to be super
# extra explicit.


positions = [walker[i][walker[i]!=0] for i in range(len(walker))]
steps = [jumps[i][jumps[i]!=0] for i in range(len(jumps))]
# Don't have to renormalized net_flux, because extra zeroes are okay.
fluxes = [sum(net_flux[i]) for i in range(len(net_flux))]
print('Overall net flux = {}'.format(sum(fluxes)))

# The following line is wrong because the arrays are the same size
# even if either array is not entirely filled.

print('Ratio of steps on energy landscape A to B = {}'.
      format(float(len(hstack(positions[0::2]))) /
             float(len(hstack(positions[1::2])))))
        


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
        counts, edges = histogram(hstack(positions[i::2]), range=(min(q), max(q)),
                                  bins=len(q), density=True)
        mids = (edges[1:] + edges[:-1]) / 2.
        ax.bar(mids, counts, color=c, edgecolor='none',
               width=mids[1] - mids[0], label='{} steps'.format(len(hstack(positions[i::2]))))
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



###########################################
# FIX
# 1. steps_executed counter
# 
#  [walker[i][walker[i]!=0] for i in range(len(walker))]


# 2. Overall netflux decomposed in graphs
# 3. Graph PDFs
###########################################