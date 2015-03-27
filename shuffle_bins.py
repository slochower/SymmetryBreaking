#!/usr/bin/python

from math import *
from numpy import *
from scipy.stats import norm
import matplotlib.cm as cm
# from scipy import optimize

execfile('../common/figures.py')

debug = True
animation, density, annotate = False, True, False


def landscapes(x, shift=0):
    return float(cos(2*x + shift * pi)*10 + 0.0)


def landscapes_noise(x, shift=0, noise=0):
    return float((cos(2*x + shift * pi)/4 + random.normal(0, 0.05 * pi, 1))/4)


def find_nearest(array, value):
    idx = (abs(array-value)).argmin()
    return array[idx]

def force(x, dx):
    return float(-(landscapes(x+dx, shift=shift) -
                   landscapes(x, shift=shift)) / dx)


#################################
# PARAMETERS
#################################
num_landscapes = 1
dx, dt, timesteps = 0.01, 1, 10000
q = arange(0, 2 * pi, dx)
shift = pi                # Shift of landscape between holo and apo

##################################
# BUILD ENERGY LANDSCAPES
##################################
E = []    # Container for energy landscapes
s = []    # List of beginning points for moving on the energy landscapes

for i in range(num_landscapes):
    if i % 2 == 0:
        j = 0.0
        j = shift
    else:
        j = shift
    E.append([landscapes(i, shift=j) for i in q])

##################################
# BROWNIAN MOTION ON THE LANDSCAPE
##################################

for i in range(num_landscapes):
    positions = []
    jumps = []
    fluxes = []
    x = pi     # Start Brownian walker at position x = pi
#    force = -diff(E[i]) / dx
    D = 0.1      # Arbitrary
#################################
    kT = 100.0  # Arbitrary
#################################
    positions.append(x)
    for t in range(timesteps):
        swap = 0
        gamma = random.normal(loc=0.0, scale=2*D*dt)
        new_x = x + (D/kT) * force(x,dx) * dt + gamma
#        new_x = x + force(x,dx) * dt
        if debug:
            print 't = {} \t x = {}, force * dt = {}, new_x = {}'.format(float(t), float(x), float(force(x,dx)*dt), float(new_x))

        if new_x > max(q):
            print 'Wrapping {} to {}'.format(new_x, min(q)+(new_x-max(q)))
            tmp = new_x - x
            new_x = min(q) + (new_x-max(q))
            swap = 1
            if debug:
                print 'Setting swap to {}'.format(swap)
                print 'Jump = {}'.format(new_x - x + swap*(max(q)-min(q)))
        elif new_x < min(q):
            if debug:
                print 'Wrapping {} to {}'.format(new_x, max(q) -
                                                 abs((new_x-min(q))))
            tmp = new_x - x
            new_x = max(q) - abs(new_x - min(q))
            swap = -1
        else:
            swap = 0
        fluxes.append(swap)
        positions.append(new_x)
        jump = new_x - x + swap*(max(q)-min(q))
        if (abs(jump) > 5):
            print new_x
            print x
            print swap
            print 'Delta = {}'.format(tmp)
            print 'Jump =  {} '.format(new_x - x + swap*(max(q)-min(q)))
            break
        jumps.append(jump)
        x = new_x
    boltzmann = exp([-E[i][j]/kT for j in range(len(E[i]))])


##################################
# PLOTTING
##################################

c = 'b'
fig, gs, axes = generate_axes_pad(nrows=3, ncols=2, v_pad=0.3,
                                  h_pad=0.2, figsize=(12, 12))
ax = axes[0][0]
ax.plot(q, E[0], color=c, lw=2)
ax.set_title('Energy landscape')

ax = axes[0][1]
ax.plot(q, [force(x, dx) for x in q], color=c, lw=2)
ax.set_title('Force')

ax = axes[1][0]
ax.plot(q, boltzmann, color='k', lw=2, label='kT = {}'.format(kT))
ax.set_title('Boltzmann distribution')
ax.legend()

ax = axes[1][1]
ax.set_title('Position density along energy landscape')
ax.plot(q, E[0], color=c, lw=2)
energies = [landscapes(j,shift=shift) for j in [positions[k] for k in
            range(len(positions))]]
ax.scatter(positions, energies,
           c='k',
           s=200, lw=0, alpha=0.003)
labels = ['Point {0}'.format(i) for i in range(len(positions))]
'''
for label, x, y in zip(labels, positions, energies):
    plt.annotate(
                 label,
                 xy = (x, y), xytext = (-20, 20),
                 textcoords = 'offset points', ha = 'right', 
                 va = 'bottom',
                 bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', 
                 alpha = 0.5),
                 arrowprops = dict(arrowstyle = '->', 
                 connectionstyle = 'arc3,rad=0'))
ax.legend()
'''
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
ax.set_title('Step sizes: $\mu = {0:.4f}, \sigma = {0:.3f}$'.format(mu,std))
ax.legend()
plt.show()
# plt.savefig('shuffle_bins_kT.{}.png'.format(float(kT)))
plt.close(fig)

if animation: 
    for k in range(timesteps):
        fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                          h_pad=0.2, figsize=(12, 12))
        ax = axes[0][0]
        cmap = cm.jet
        c = linspace(min(positions), max(positions), len(positions))
        ax.set_title('Energy landscape')
        ax.plot(q, E[0], color=c, lw=2)
        cax = ax.scatter(positions[k],
                         E[0][where(q == find_nearest(q, positions[k]))[0][0]],
                         #c=c[k], cmap=cmap,
                         c='k',
                         s=200, lw=0, alpha=0.5,
                         label='Net boundary crossings = {} \n Total boundary'
                         ' crossings = {}'
                         .format(sum(fluxes), sum(absolute(fluxes))))

        ax.legend()
        plt.savefig('tmp.{:02d}.png'.format(k))
        plt.close(fig)
    from subprocess import call
    call(["ffmpeg", "-framerate", "90", "-i", "tmp.%02d.png", "-c:v",
          "libx264", "-r", "30", "-pix_fmt", "yuv420p", "motion.mp4"])
    call(["rm", "tmp*.png"])
    import glob
    import os
#    for fl in glob.glob("tmp*.png"):
#        os.remove(fl)

if density: 
    fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))
    ax = axes[0][0]
    ax.set_title('Energy landscape')
    ax.plot(q, E[0], color=c, lw=2)
    ax.scatter(positions, [-15 for k in range(len(positions))],
                         c='k',
                         s=200, lw=0, alpha=0.005,
                         label='Net boundary crossings = {} \n Total boundary'
                         ' crossings = {}'
                         .format(sum(fluxes), sum(absolute(fluxes))))
    ax.legend()
    plt.savefig('shuffle_bins_density.png')
    plt.close(fig)

if annotate:
    fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                      h_pad=0.2, figsize=(12, 12))
    ax = axes[0][0]
    ax.set_title('Energy landscape')
    ax.plot(q, E[0], color=c, lw=2)
    energies = [landscapes(j) for j in [positions[k] for k in
                range(len(positions))]]
    ax.scatter(positions, energies,
               c='k',
               s=200, lw=0, alpha=1)
    labels = ['Point {0}'.format(i) for i in range(len(positions))]
    for label, x, y in zip(labels, positions, energies):
        plt.annotate(
                     label,
                     xy = (x, y), xytext = (-20, 20),
                     textcoords = 'offset points', ha = 'right', 
                     va = 'bottom',
                     bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', 
                     alpha = 0.5),
                     arrowprops = dict(arrowstyle = '->', 
                     connectionstyle = 'arc3,rad=0'))
    ax.legend()
    plt.show()
    plt.close(fig)
