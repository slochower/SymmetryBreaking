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

    positions = []
    fluxes = []
    # Crossings of the flux barrier
    # -1: left crossing of the boundary
    #  0: no crossing this step
    # +1: right crossing of the boundary
    flux_point = pi

    if (i % 2) == 0: 
        state = 0
    else:
        state = 1

    if debug:
        print('\nStarting walker at x = {}'.format(x))
        print('Running for a maximum of {} steps on this landscape'.format(steps_on_this_landscape))
        print('Recording position = {}'.format(x))
    # Record the initial position of the walker at t = 0 and set
    # the flux at this time to be 0, by definition.
    positions.append(x)
    fluxes.append(0)
    # Since we have recorded the position and flux, count that as a 
    # time step.
    t = 1
    # Each iteration through the loop adds two timesteps: one for MC and
    # one for BD, so if there is only 1 step remaining, then quit early.
    quit_early = False
    if steps_on_this_landscape == 1:
        quit_early = True
    while t < steps_on_this_landscape-1:
        #######################################
        # BROWNIAN WALK
        #######################################
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(q[:-1], forces[state], x)
        new_x = x + (D / kT) * F * dt + g

        if (abs(new_x - x) > max(q)):
            raise Exception('Jumps are too big.')

        #######################################
        # MEASURE FLUX
        #######################################
        if (x < flux_point) and (new_x >= flux_point):
            swap = 1
        elif (x > flux_point) and (new_x <= flux_point):
            swap = -1
        else:
            swap = 0

        ########################################
        # KEEP TRACK OF PBC CROSSINGS
        ########################################
        if new_x > max(q):
            new_x = min(q) + (new_x - max(q))
        elif new_x < min(q):
            new_x = max(q) - abs(new_x - min(q))
        if debug:
            print('t = {}, landscape = {}, x = {}, new_x = {}, status = {}'.
                  format(t, state, x, new_x, str('BD')))
            print('Recording position = {}'.format(new_x))

        ####################
        # RECORD KEEPING
        ####################
        t += 1
        fluxes.append(swap)
        positions.append(new_x)
        x = new_x
        if quit_early: break
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
        if (MC is False):
            raise Execption('Brownian dynamics without MC is not implemented.')

        if (MC is True and ((t+1) % MC_interval) == 0):
            E_transition = energy_lookup(q, energies[(state+1) % 2], x)
            E_now = energy_lookup(q, energies[state], x)
            delta = E_transition - E_now
            # If delta < 0, step to next landscape
            if (delta < 0):
                if debug:
                    print('t = {}, landscape = {}, x = {}, status = {}'.
                          format(t, state, x, str('MC ACCEPT (delta E)')))
                #    print('Recording position = {}'.format(x))
                #positions.append(x)
                #fluxes.append(0)
                #t += 1
                break
            # If delta !< 0, compute p_accept and pick a random number.
            p_accept = exp(-delta/kT)
            r = random.random()
            if (p_accept > r):
                if debug:
                    print('t = {}, landscape = {}, x = {}, status = {}'.
                          format(t, state, x, str('MC ACCEPT (p_accept > r)')))
                #    print('Recording position = {}'.format(x))
                #positions.append(x)
                #fluxes.append(0)
                #t += 1
                break
            # If p_accept !> r, append the current position and 
            # then take a Brownian dynamics step.
            if debug:
                print('t = {}, landscape = {}, x = {}, status = {}'.
                      format(t, state, x, str('MC FAIL')))
                print('Recording position = {}'.format(x))
            t += 1
            positions.append(x)
            fluxes.append(0)
            pass

    return(positions, fluxes)


#################################
# PARAMETERS
#################################
debug = False
MC = True
flashing = False
plot_together = True
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
MC_interval = 1
# Number of steps between Monte Carlo iterations
total_timesteps = 100000
# Cumulative number of steps   for the simulation (across landscapes)
D = 0.01
# Arbitrary -- these two work together!
kT = 10
# Arbitrary -- these two work together!
shift = [0, pi]


###################################
# PRE-CALCULATE ENERGY, FORCE, PDF
###################################
energies = [energy(q, shift[i]) for i in range(len(shift))]
if flashing:
    energies[1] = array([mean(energies[0])]*len(energies[0]))
forces = [[force(k, dx, energies[i]) for k in range(len(q)-1)]
          for i in range(len(energies))]
boltzmann = [exp(-energies[i]/kT) for i in range(len(energies))]
pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(len(boltzmann))]

##################################
# SIMULATE
##################################
# N.B. total_timesteps is the maximum length that either dimension could be..
# Infinitessimal chance the walker will be exactly zero for the last step,
# so we can simply trim zeros from the array to condense the data later.
walker, net_flux = zeros((2, total_timesteps)), zeros((2, total_timesteps))
steps_executed, landscapes_sampled, start = 0, 0, 0.0
# Being really really explicit right now is necessary for understanding...
steps_on_A, steps_on_B = 0, 0
start_timer = timer()

while steps_executed < total_timesteps:
    this_run = simulate(start, dx, D, kT, dt, shift, forces, energies,
                        landscapes_sampled, 
                        min(total_timesteps, total_timesteps - steps_executed))
    
    on_A = (landscapes_sampled+1) % 2
    if on_A:
        walker[0][steps_on_A:steps_on_A+len(this_run[0])] = this_run[0]
        net_flux[0][steps_on_A:steps_on_A+len(this_run[1])] = this_run[1]
        steps_on_A += len(this_run[0])
        start = walker[0][steps_on_A-1]
    else:
        walker[1][steps_on_B:steps_on_B+len(this_run[0])] = this_run[0]
        net_flux[1][steps_on_B:steps_on_B+len(this_run[1])] = this_run[1]
        steps_on_B += len(this_run[0])
        start = walker[1][steps_on_B-1]
    
    steps_executed += len(this_run[0])
    landscapes_sampled += 1
    
    if debug:
        print('Total steps = {}'.format(steps_executed))

    
simulation_time = timer() - start_timer
print('###############################################################')
print('Simulation took {} seconds'.format(simulation_time))
print('Total landscapes sampled: {}'.format(landscapes_sampled))

################################
# CONDENSE THE DATA
################################
A = trim_zeros(walker[0], 'b')
B = trim_zeros(walker[1], 'b')
A_flux = trim_zeros(net_flux[0], 'b')
B_flux = trim_zeros(net_flux[1], 'b')

print('Overall net flux = {}'.format(sum(A_flux)+sum(B_flux)))
print('Ratio of steps on energy landscape A to B = {}'.
      format(float(len(A)) / float(len(B))))
print('###############################################################')

if plot_together:
    fig, gs, axes = generate_axes_pad(nrows=2, ncols=2, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))

    c = clrs[0]
    ax = axes[0][0]
    counts, edges = histogram(A, range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    ax.bar(mids, counts, color=c, edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(A)))
    ax.plot(q, pdf[0], color='k', lw=2)
    ax.legend()

    ax = axes[1][0]
    ax.step(range(len(A_flux)), cumsum(A_flux), color=c, lw=2)
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Steps on this landscape')
    ax.set_title('Net flux this landscape = {}\nTotal crossings = {}'.format(sum(A_flux), sum(A_flux != 0)))

    c = clrs[1]
    ax = axes[0][1]
    counts, edges = histogram(B, range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    ax.bar(mids, counts, color=c, edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(B)))
    ax.plot(q, pdf[1], color='k', lw=2)
    ax.legend()

    ax = axes[1][1]
    ax.step(range(len(B_flux)), cumsum(B_flux), color=c, lw=2)
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Steps on this landscape')
    ax.set_title('Net flux this landscape = {}\nTotal crossings = {}'.format(sum(B_flux), sum(B_flux != 0)))

    
    fig.text (0.5, 0.98, 'Total landscapes = {}, MC interval = {}\nInitial Position = {}'.format(landscapes_sampled, MC_interval, float(0)), horizontalalignment='center',
         verticalalignment='top', size=fs['title'])
    plt.show()

if interpolation_tests:
    execfile('2D_walk.py')

if parameter_scan is True:
    execfile('parameter_scan.py')
