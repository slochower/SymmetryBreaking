#!/usr/bin/python

from math import *
from numpy import *
from scipy import ndimage
from timeit import default_timer as timer
import time

execfile('../common/figures.py')

routine = ['BDMC', 'MCMC', 'BD2D'][-1]

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




#################################
# PARAMETERS
#################################
debug, debug_alternative, log = False, False, False
MC = True
flashing = False
plot_together = True

if log == True:
    import sys
    old_stdout = sys.stdout
    log_file = open('shuffle_bins_log', 'w')
    sys.stdout = log_file

dx = 0.01
# Reaction coordinate spacing (necessary for interpolation)
q = arange(0, 2 * pi, dx)
# Define reaction coordinate
dt = 1
# Time resolution (do not change this)
D = 0.01
# Arbitrary -- these two work together!
kT = 10
# Arbitrary -- these two work together!
shift = [0, pi]
# Offset between energy landscapes
MC_interval = 1
# Number of steps between Monte Carlo iterations (for moving landscapes)
total_timesteps = 1000000
# Cumulative number of steps for the simulation (across landscapes)

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
if 'BDMC' in routine:
    print('Running mixed BD/MC')
    execfile('BDMC.py')
if 'MCMC' in routine:
    print('Running MC/MC')
    execfile('MCMC.py')

if 'BD2D' in routine:
    print('Running 2D BD')
    execfile('2D_walk.py')



while steps_executed < total_timesteps:
    if 'BDMC' in routine:
        this_run = simulateBDMC(start, dx, D, kT, dt, shift, forces, energies,
                        landscapes_sampled, 
                        min(total_timesteps, total_timesteps - steps_executed))
    if 'MCMC' in routine:
        this_run = simulateMCMC(start, dx, D, kT, dt, shift, forces, energies,
                        landscapes_sampled, 
                        min(total_timesteps, total_timesteps - steps_executed))
        # print ('landscapes sampled = {}'.format(landscapes_sampled))
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
if len(B) > 0:
    print('Ratio of steps on energy landscape A to B = {}'.
          format(float(len(A)) / float(len(B))))
print('###############################################################')
if log == True:
    sys.stdout = old_stdout
    log_file.close()

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
    ax.step([float(i)/1000 for i in range(len(A_flux))], cumsum(A_flux), color=c, lw=2)
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Steps on this landscape (thousands)')
    ax.grid()
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
    ax.step([float(i)/1000 for i in range(len(B_flux))], cumsum(B_flux), color=c, lw=2)
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Steps on this landscape (thousands)')
    ax.grid()
    ax.set_title('Net flux this landscape = {}\nTotal crossings = {}'.format(sum(B_flux), sum(B_flux != 0)))

    
    fig.text (0.5, 0.98, 'Total landscapes = {}, MC interval = {}\nInitial Position = {}'.format(landscapes_sampled, MC_interval, float(0)), horizontalalignment='center',
         verticalalignment='top', size=fs['title'])
    plt.show()

