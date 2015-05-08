#!/usr/bin/python

from math import *
from numpy import *
from scipy import ndimage
from timeit import default_timer as timer
import time


routine = ['BDMC', 'MCMC', 'BD2D'][-1]
if 'BDMC' in routine:
    options = ['2D_MC']
if 'BD2D' in routine:
    options = []
plot = True
debug, debug_alternative, log = False, False, False


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
if log == True:
    import sys
    old_stdout = sys.stdout
    log_file = open('shuffle_bins_log', 'w')
    sys.stdout = log_file

if 'BDMC' in routine:
    print('Running mixed BD/MC')
    # Only implemented for two 1D landscapes.
    dx = 0.01
    q = arange(0, 2 * pi, dx)
    dt = 1
    D = 0.01
    kT = 10
    shift = [0, pi]
    MC_interval = 1
    total_timesteps = 1000000
    MC = True
    flashing = False

    execfile('BDMC.py')
    if '2D_MC' in options:
        print('MC allowed to work on either surface.')
        MC_side_attempts = 0
        MC_orthogonal_attempts = 0

    initialize_BDMC()
    start_timer = timer()
    while steps_executed < total_timesteps:
        if '2D_MC' in options:
            this_run = simulate_BDMC_both_surfaces(start, dx, D, kT, dt, shift, forces, energies, landscapes_sampled, min(total_timesteps, total_timesteps - steps_executed))
            MC_orthogonal_attempts += this_run[2]
            MC_side_attempts += this_run[3]
        else:
            this_run = simulate_BDMC(start, dx, D, kT, dt, shift, forces,
                                    energies, landscapes_sampled,
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
    simulation_time = timer() - start_timer
    print('###############################################################')
    print('Simulation took {} seconds'.format(simulation_time))
    print('Total landscapes sampled: {}'.format(landscapes_sampled))

    A = trim_zeros(walker[0], 'b')
    B = trim_zeros(walker[1], 'b')
    A_flux = trim_zeros(net_flux[0], 'b')
    B_flux = trim_zeros(net_flux[1], 'b')

    print('Overall net flux = {}'.format(sum(A_flux)+sum(B_flux)))
    if len(B) > 0:
        print('Ratio of steps on energy landscape A to B = {}'.
              format(float(len(A)) / float(len(B))))
    print('###############################################################')

if 'MCMC' in routine:
    print('Running MC/MC')
    # Only implemented for two 1D landscapes.
    dx = 0.01
    q = arange(0, 2 * pi, dx)
    dt = 1
    D = 0.01
    kT = 10
    shift = [0, pi]
    MC_interval = 1
    total_timesteps = 1000000
    MC = True
    flashing = False

    execfile('MCMC.py')

    energies, forces, boltzmann, pdf, walker, net_flux, steps_executed, landscapes_sampled, start, steps_on_A, steps_on_B = initialize_MCMC()
    start_timer = timer()
    while steps_executed < total_timesteps:
        this_run = simulate_MCMC(start, dx, D, kT, dt, shift, forces,
                                energies, landscapes_sampled,
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
    simulation_time = timer() - start_timer
    print('###############################################################')
    print('Simulation took {} seconds'.format(simulation_time))
    print('Total landscapes sampled: {}'.format(landscapes_sampled))

    A = trim_zeros(walker[0], 'b')
    B = trim_zeros(walker[1], 'b')
    A_flux = trim_zeros(net_flux[0], 'b')
    B_flux = trim_zeros(net_flux[1], 'b')

    print('Overall net flux = {}'.format(sum(A_flux)+sum(B_flux)))
    if len(B) > 0:
        print('Ratio of steps on energy landscape A to B = {}'.
              format(float(len(A)) / float(len(B))))
    print('###############################################################')
    
if 'BD2D' in routine:
    print('Running 2D BD')
    dx, dy = 0.01, 0.01
    q = arange(0, 2 * pi, dx)
    transition_distance = 40.0
    y_one = [0]*len(q)
    y_two = [0+transition_distance]*len(q)
    r = mgrid[y_one[-1]:y_two[-1]:100j]
    dt = 1
    D = 0.01
    kT = 10
    total_timesteps = 10

    MC = False
    flashing = False

    execfile('2D_walk.py')
    energy, force, boltzmann, pdf, walker, net_flux, steps_executed, landscapes_sampled, start, steps_on_A, steps_on_B = initialize_2DBD()
    start_timer = timer()
    while steps_executed < total_timesteps:
        this_run = simulate_2DBD(start[0], dx, start[1], dy, D, kT, dt, force,
                                energy, min(total_timesteps, total_timesteps - steps_executed))
        walker[steps_on_A:steps_on_A+len(this_run[0])] = this_run[0]
        net_flux[steps_on_A:steps_on_A+len(this_run[1])] = this_run[1]
        steps_on_A += len(this_run[0])

        
        steps_executed += len(this_run[0])
        landscapes_sampled += 1        
    simulation_time = timer() - start_timer
    print('###############################################################')
    print('Simulation took {} seconds'.format(simulation_time))
    print('Total landscapes sampled: {}'.format(landscapes_sampled))

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

if plot:
    execfile('../common/figures.py')

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
