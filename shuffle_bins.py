#!/usr/bin/python

from math import *
from numpy import *
from scipy import ndimage
from timeit import default_timer as timer
import time

execfile('../common/figures.py')

routine = ['BDMC']

def energy(x, shift=0):
    period = 2 * pi
    amplitude = 40
    window = 20
    sawtooth = amplitude * ((x + shift) / period -
                            floor(1 / 2 + (x + shift) / period))

    padded = concatenate((sawtooth, sawtooth, sawtooth))
    smoothed = gaussian_filter(padded, window)
    return(smoothed[len(sawtooth):2*len(sawtooth)])


def gaussian_filter(x, N):
    return(ndimage.filters.gaussian_filter1d(input=x, sigma=N))


def force(i, dx, array):
    F = -((array[i+1] - array[i]) / dx)
    return(F)


def force_lookup(q, array, value):
    push = interp(value, q, array)
    return(push)

def energy_lookup(position, array, value):
    light = interp(value, position, array)
    return(light)

# @profile
def simulate(x, dx, D, kT, dt, shift, forces, energies, i, steps_on_this_landscape, **kwargs):
    ##################################
    # BROWNIAN MOTION ON A LANDSCAPE
    ##################################
    positions = []
    jumps = []
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
    flux_point = pi

    #################
    # This could be carried out before calling simulate()
    if (i % 2) == 0: 
        state = 0
    else:
        state = 1
    ################
    if debug:
        print('Starting walker at x = {}'.format(x))
        print('Running for a maximum of {} steps on this landscape'.format(steps_on_this_landscape))

    for t in range(steps_on_this_landscape):
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
        if (MC is False):
            print('This behavior might be incorrect because the code'\
                  ' will still run steps_on_this_landscape and then switch.')

        if (MC is True and (t % MC_interval) == 0 and t > MC_interval):
            if debug:
                print('MC check')
            E_transition = energy_lookup(q, energies[(state+1) % 2], x)
            E_now = energy_lookup(q, energies[state], x)
            delta = E_transition - E_now
            p_accept = exp(-delta/kT)
            r = random.random()
            if (delta < 0):
                if debug:
                    print('Move because delta < 0 after {} steps (not including MC step).'.format(t))
                reason = 0
                break
            if (p_accept > r):
                if debug:
                    print('Move because p_accept > r after {} steps (not including MC step)'.format(t))
                reason = 1
                break

        #######################################
        # BROWNIAN WALK
        #######################################
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(q[:-1], forces[state], x)
        new_x = x + (D / kT) * F * dt + g
        ########################################
        # KEEP TRACK OF PBC CROSSINGS
        ########################################
        if new_x > max(q):
            new_x = min(q) + (new_x - max(q))
            wrap = 1
        elif new_x < min(q):
            new_x = max(q) - abs(new_x - min(q))
            wrap = -1
        else:
            wrap = 0
        if debug:
            print 'landscape = {} \t t = {} \t x = {}, force * dt = {}, new_x = {}' \
                  .format(state, float(t), float(x),
                          float(F * dt), float(new_x))
        #######################################
        # MEASURE FLUX
        #######################################
        if (x < flux_point) and (new_x >= flux_point):
            swap = 1
        elif (x > flux_point) and (new_x <= flux_point):
            swap = -1
        else:
            swap = 0

        ####################
        # CLEANUP AND RECORD
        ####################
        fluxes.append(swap)
        positions.append(new_x)

        jump = new_x - x + wrap * (max(q) - min(q))
        if (abs(jump) > max(q)):
            raise Exception('Jumps are too big.')
        jumps.append(jump)
        x = new_x

    if (t == (timesteps - 1) and MC is True):
        raise Exception('There were no transitions between landscapes'\
                        ' during this simulation.'\
                        ' We are deviating from what we know and'\
                        ' understand.')

    if plot_every is True:
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        ax.plot(q, energies[state], color='k', lw=2)
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
                             [energy_lookup(q, energies[state], positions[k])
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
    # Returning t prevents recording MC step as a step, as MC step will become
    # initial position on next landscape.
    # print('Actually ran for {} steps'.format(t+1))
    return(positions, jumps, fluxes, t)


#################################
# PARAMETERS
#################################
debug, deterministic, MC = False, False, True
flashing = False
plot_together, plot_every = True, False
interpolation_tests, parameter_scan = False, False

dx = 0.01
# Reaction coordinate spacing, although the walker is basically 
# a continuous variable.
q = arange(0, 2 * pi, dx)
# Define reaction coordinate
dt = 1
# Time resolution
# Don't change this because it is tied in with the array spacing
# e.g. positions[1] - positions[0] is the step from time 0 to time 1.
max_MC_attempts = 100
# Need a guess for MC acceptance rate to pre-allocate arrays
MC_interval = 1
# Number of steps between Monte Carlo iterations
timesteps = max_MC_attempts*MC_interval
# Above, maximum the number of time steps to sample per landscape
total_timesteps = 10
# Cumulative number of steps   for the simulation (across landscapes)
D = 0.01
# Arbitrary -- these two work together!
kT = 10
# Arbitrary -- these two work together!
shift = [0, pi]


#################################
# SIMULATE
#################################


energies = [energy(q, shift[i]) for i in range(len(shift))]
if flashing:
    energies[1] = array([mean(energies[0])]*len(energies[0]))
forces = [[force(k, dx, energies[i]) for k in range(len(q)-1)]
          for i in range(len(energies))]
boltzmann = [exp(-energies[i]/kT) for i in range(len(energies))]
pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(len(boltzmann))]

walker = empty((total_timesteps/MC_interval, timesteps))
net_flux = empty((total_timesteps/MC_interval, timesteps))
steps_executed = []
landscapes_sampled = 0
start = timer()
while sum(steps_executed) < total_timesteps:
    # print 'Landscape {}'.format(landscapes_sampled)
    timesteps_remaining = int(total_timesteps - sum(steps_executed))
    if sum(steps_executed) % 10000 == 0:
        print '{} steps remaining'.format(timesteps_remaining)
    if timesteps_remaining > timesteps:
        steps_on_this_landscape = timesteps
    else:
        steps_on_this_landscape = timesteps_remaining
    if landscapes_sampled == 0:
        position = 0
    else:
        position = walker[landscapes_sampled-1][steps_executed[-1]-2]
    this_run = simulate(position, dx, D, kT, dt, shift, forces, energies, landscapes_sampled, steps_on_this_landscape)
    walker[landscapes_sampled,0:len(this_run[0])] = this_run[0]
    net_flux[landscapes_sampled,0:len(this_run[2])] = this_run[2]
    steps_executed.append(this_run[3]+1)
    # print('Recorded running for {} steps'.format(steps_executed[-1]))
    landscapes_sampled += 1
    tmp = [walker[i][0:steps_executed[i]-1] for i in range(landscapes_sampled)]
    counts, edges = histogram(hstack(tmp[0::2]), range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    plt.figure()
    plt.bar(mids, counts, color=clrs[0], edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(hstack(tmp[0::2]))))
    plt.plot(q, pdf[0], color='k', lw=2)
    plt.ylim([0,1])
    plt.savefig('landscape_{}.png'.format(str(landscapes_sampled)))
    plt.close()

simulation_time = timer() - start
total_landscapes = landscapes_sampled
print('Simulation took {} seconds'.format(simulation_time))
print('Total landscapes sampled: {}'.format(total_landscapes))

################################
# CONDENSE THE DATA
################################shape
positions = [walker[i][0:steps_executed[i]-1] for i in range(total_landscapes)]
# Hack to save the final point
positions[total_landscapes-1] = hstack([positions[total_landscapes-1],  walker[total_landscapes-1][steps_executed[-1]-1]])
fluxes = [net_flux[i][0:steps_executed[i]-1] for i in range(total_landscapes)]
# Hack to save the final point
fluxes[total_landscapes-1] = hstack([fluxes[total_landscapes-1],  net_flux[total_landscapes-1][steps_executed[-1]-1]])

print('Overall net flux = {}'.format(sum(hstack(fluxes))))
print('Ratio of steps on energy landscape A to B = {}'.
      format(float(len(hstack(positions[0::2]))) /
             float(len(hstack(positions[1::2])))))
print('Monte Carlo acceptance ratio is approximately 1 out of every {} attempts'.
      format((mean(steps_executed)-2)/MC_interval))
# Minus 2 because t0 --> t1 doesn't have an MC step and
# the final step is an MC step. This will give MC attempts/MC interval.
print('###############################################################')

##########################
# CLEANUP
##########################
del walker
del net_flux
# Needed for plotting:
# del positions
# del fluxes

if plot_together:
    fig, gs, axes = generate_axes_pad(nrows=2, ncols=2, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))

    c = clrs[8]
    ax = axes[0][0]
    counts, edges = histogram(hstack(positions[0::2]), range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    ax.bar(mids, counts, color=c, edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(hstack(positions[0::2]))))
    ax.plot(q, pdf[0], color='k', lw=2)
    ax.legend()

    ax = axes[1][0]
    ax.step(range(len(hstack(fluxes[0::2]))), cumsum(hstack(fluxes[0::2])), color=c, lw=2)
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Steps on this landscape')
    ax.set_title('Net flux this landscape = {}\nTotal crossings = {}'.format(sum(hstack(fluxes[0::2])), sum((hstack(fluxes[0::2])) != 0)))
    
    ax = axes[0][1]
    counts, edges = histogram(hstack(positions[1::2]), range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    ax.bar(mids, counts, color=c, edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(hstack(positions[1::2]))))
    ax.plot(q, pdf[1], color='k', lw=2)
    ax.legend()

    ax = axes[1][1]
    ax.step(range(len(hstack(fluxes[1::2]))), cumsum(hstack(fluxes[1::2])), color=c, lw=2)
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Steps on this landscape')
    ax.set_title('Net flux this landscape = {}\nTotal crossings = {}'.format(sum(hstack(fluxes[1::2])), sum((hstack(fluxes[1::2])) != 0)))

    
    fig.text (0.5, 0.98, 'Total landscapes = {}, MC interval = {}\nInitial Position = {}'.format(total_landscapes, MC_interval, float(0)), horizontalalignment='center',
         verticalalignment='top', size=fs['title'])
    plt.show()

if interpolation_tests:
    execfile('2D_walk.py')

if parameter_scan is True:
    execfile('parameter_scan.py')
