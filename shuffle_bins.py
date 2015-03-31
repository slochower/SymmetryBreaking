#!/usr/bin/python

from math import *
from numpy import *
from scipy.stats import norm
# from scipy import optimize
import time
import datetime

execfile('../common/figures.py')

debug, deterministic = False, False                 # Tests
animation, density, annotate = False, True, False   # Plotting options
plot_separate = False                               # Plotting options
plot_together = True


def landscapes(x, shift=0):
    # return float(cos(2*x + shift * pi)*10 + 0.0)
    period = 1.8*pi
    amplitude = 10
    return(float(amplitude * ((x+shift)/period -
                 floor(1/2 + (x+shift)/period))))


def landscapes_noise(x, shift=0, noise=0):
    return float((cos(2*x + shift * pi)/4 + random.normal(0, 0.05 * pi, 1))/4)


def find_nearest(array, value):
    idx = (abs(array-value)).argmin()
    return array[idx]


def force(x, dx, s):
    return float(-(landscapes(x+dx, shift=s) -
                   landscapes(x, shift=s)) / dx)


#################################
# PARAMETERS
#################################
num_landscapes = 5          # Number of protein states
dx = 0.01                   # Reaction coordinate spacing
q = arange(0, 2 * pi, dx)   # Define reaction coordinate
dt = 1                      # Time resolution (not really used)
timesteps = 100000           # Number of steps for walker
E = []                      # Container for energy landscapes
walker = []                 # Positions of the walker over all possible landscapes
net_flux = []

if plot_together:
    fig, gs, axes = generate_axes_pad(nrows=5, ncols=3, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))


for i in range(num_landscapes):
    if i % 2 == 0:
        shift = 0.0
        c = clrs[i % 9]
    else:
        shift = pi
        c = clrs[i % 9]
    E.append([landscapes(pos, shift=shift) for pos in q])

    ##################################
    # BROWNIAN MOTION ON THE LANDSCAPE
    ##################################

    positions = []        # List of walker positions
    jumps = []            # List of walker steps
    fluxes = []           # List of PBC transitions (+/- 1)
    if i == 0:        # Start Brownian walker at position pi
        x = pi
    D = 0.1      # Arbitrary -- these two work together!
    kT = 100.0   # Arbitrary -- these two work together!
    positions.append(x)   # Set initial position
    walker.append(x)
    print('Starting walker at x = {}'.format(x))

    if debug:
        print('Starting walker at x = {}'.format(x))
    for t in range(timesteps):
        swap = 0
        gamma = random.normal(loc=0.0, scale=2*D*dt)
        new_x = x + (D/kT) * force(x, dx, shift) * dt + gamma
        if deterministic:
            new_x = x + force(x, dx, shift) * dt
        if debug:
            print 't = {} \t x = {}, force * dt = {}, new_x = {}' \
                  .format(float(t), float(x),
                          float(force(x, dx, shift)*dt), float(new_x))

        if new_x > max(q):
            new_x = min(q) + (new_x-max(q))
            swap = 1
            if debug:
                print 'Wrapping {} to {}' \
                      .format(new_x, min(q)+(new_x-max(q)))
                print 'Setting swap to {}'.format(swap)
                print 'Jump = {}'.format(new_x - x + swap*(max(q)-min(q)))
        elif new_x < min(q):
            if debug:
                print 'Wrapping {} to {}'.format(new_x, max(q) -
                                                 abs((new_x-min(q))))
            new_x = max(q) - abs(new_x - min(q))
            swap = -1
        else:
            swap = 0
        fluxes.append(swap)
        positions.append(new_x)
        walker.append(new_x)
        jump = new_x - x + swap*(max(q)-min(q))
        if (abs(jump) > max(q)):
            print new_x
            print x
            print swap
            print 'Jump =  {} '.format(new_x - x + swap*(max(q)-min(q)))
            print 'WARNING: JUMPED ACROSS AN ENTIRE LANDSCAPE.'
        jumps.append(jump)
        x = new_x
    boltzmann = exp([-E[i][j]/kT for j in range(len(E[i]))])
    pdf = boltzmann / sum([boltzmann[j]*dx for j in boltzmann])

    print 'Net flux on energy landscape {} = {} in {} steps' \
          .format(i, sum(fluxes), timesteps)
    net_flux.append(fluxes)
    print 'Final x = {}'.format(x)

    if plot_separate:

        ##################################
        # PLOTTING INDIVIDUAL LANDSCAPES
        ##################################

        fig, gs, axes = generate_axes_pad(nrows=3, ncols=2, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        ax.plot(q, E[i], color=c, lw=2)
        ax.set_title('Energy landscape')

        ax = axes[0][1]
        ax.plot(q, [force(x, dx, shift) for x in q], color=c, lw=2)
        ax.set_title('Force')

        ax = axes[1][0]
        ax.plot(q, pdf, color='k', lw=2, label='kT = {}'.format(kT))
        ax.set_title('Probability density function')
        ax.legend()

        ax = axes[1][1]
        ax.set_title('Position density along energy landscape')
        ax.plot(q, E[i], color=c, lw=2)
        energies = [landscapes(j, shift=shift) for j in [positions[k] for k in
                    range(len(positions))]]
        ax.scatter(positions, energies,
                   c='k',
                   s=200, lw=0, alpha=100*(1/float(timesteps)))

        ax = axes[2][0]
        counts, edges = histogram(positions, range=(min(q), max(q)),
                                  bins=10, normed=True)
        mids = (edges[1:]+edges[:-1])/2.
        ax.bar(mids, counts, color=c, width=mids[1]-mids[0])

        ax.set_title('Histogram of positions')
        ax.legend()

        ax = axes[2][1]
        counts, edges = histogram(jumps, range=(min(jumps), max(jumps)),
                                  bins=10, normed=True)
        mids = (edges[1:]+edges[:-1])/2.
        ax.bar(mids, counts, color=c, width=mids[1]-mids[0])
        mu, std = norm.fit(jumps)
        xmin, xmax = min(jumps), max(jumps)
        x = linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        ax.plot(x, p, 'k', linewidth=2)
        ax.set_title('Step sizes: $\mu = {0:.4f}, '
                     '\sigma = {0:.3f}$'.format(mu, std))
        ax.legend()
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')
        # plt.savefig('{}-individual-{}.png'.format(st,ts)
        # plt.show()
        plt.savefig('{}-individual.png'.format(i+1))
        # plt.close(fig)


        ##################################
        # PLOTTING MULTIPLE LANDSCAPES
        ##################################
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        counts, edges = histogram(walker, range=(min(q), max(q)),
                                  bins=10, normed=False)
        mids = (edges[1:]+edges[:-1])/2.
        ax.bar(mids, counts, color=c, width=mids[1]-mids[0],
               label='Net flux = {} \n Total steps = {}'.
               format(sum(net_flux), shape(net_flux)[0]*shape(net_flux)[1]))

        ax.set_title('Histogram of Brownian walker positions, after {} landscapes'.format(i+1))
        ax.legend()
        #plt.show()
        plt.savefig('{}.png'.format(i+1))

    if plot_together:
        ##################################
        # PLOTTING LANDSCAPES TOGETHER
        ##################################
        if i < 5:

            ax = axes[i][0]
            ax.plot(q, E[i], color=c, lw=2)
            ax.scatter(positions[0], landscapes(positions[0], shift=shift),
                       color='k', s=40, zorder=10, label='Start')
            ax.scatter(positions[-1], landscapes(positions[-1], shift=shift),
                       color='k', alpha=0.5, s=40, zorder=10, label='End')
            # ax.legend()
            ax.set_title('Energy landscape')

            ax = axes[i][1]
            counts, edges = histogram(positions, range=(min(q), max(q)),
                                      bins=10, normed=True)
            mids = (edges[1:]+edges[:-1])/2.
            ax.bar(mids, counts, color=c, width=mids[1]-mids[0])
            ax.legend()
            ax.set_title('Histogram, Net flux = {}'.format(sum(fluxes)))
            # Add label

            ax = axes[i][2]
            counts, edges = histogram(walker, range=(min(q), max(q)),
                                      bins=10, normed=False)
            mids = (edges[1:]+edges[:-1])/2.
            ax.bar(mids, counts, color=c, width=mids[1]-mids[0])

            ax.set_title('Cumulative, Net flux = {}'.format(sum(net_flux)))
            ax.legend()
            plt.hold()
plt.show()


print('Overall net flux = {}'.format(sum(net_flux)))


if animation:
    for k in range(timesteps):
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        c = linspace(min(positions), max(positions), len(positions))
        ax.set_title('Energy landscape')
        ax.plot(q, E[i], color=c, lw=2)
        ax.scatter(positions[k],
                   E[i][where(q == find_nearest(q,
                        positions[k]))[0][0]],
                   c='k',
                   s=200, lw=0, alpha=0.5,
                   label='Net boundary crossings = {} \n '
                         'Total boundary'
                         ' crossings = {}'
                         .format(sum(fluxes),
                                 sum(absolute(fluxes))))

        ax.legend()
        plt.savefig('tmp.{:02d}.png'.format(k))
        plt.close(fig)
    from subprocess import call
    call(["ffmpeg", "-framerate", "90", "-i", "tmp.%02d.png", "-c:v",
          "libx264", "-r", "30", "-pix_fmt", "yuv420p", "motion.mp4"])
