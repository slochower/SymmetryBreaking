#!/usr/bin/python

from math import *
from numpy import *
import matplotlib.cm as cm
# from scipy import optimize

execfile('../common/figures.py')


def landscapes(x, shift=0, noise=0):
#    return float((cos(2*x + shift * pi)/4 + random.normal(0, 0.05 * pi, 1))/4)
    return float(cos(2*x + shift * pi)*10 + 0.0)
#    return float(random.normal(0, 0.05 * pi, 1))

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
    else:
        j = shift
    E.append([landscapes(i, shift=j) for i in q])

##################################
# BROWNIAN MOTION ON THE LANDSCAPE
##################################
# force = [-diff(E[i])/dx for i in range(num_landscapes)]
for i in range(num_landscapes):
    positions = []
    jumps = []
    fluxes = []
    x = pi     # Start Brownian walker at position x = pi
    force = -diff(E[i]) / dx
    D = 0.1      # Arbitrary
    kT = 1.0  # Arbitrary
    for t in range(timesteps):
        swap = 0
        gamma = random.normal(loc=0.0, scale=2*D*dt)
        new_x = x + (D/kT) * force[x] * dt + gamma
        if new_x > max(q):
#            print 'Wrapping {} to {}'.format(new_x, min(q)+(new_x-max(q)))
            tmp = new_x - x
            new_x = min(q) + (new_x-max(q))
            swap = 1
#            print 'Setting swap to {}'.format(swap)
#            print 'Jump = {}'.format(new_x - x + swap*(max(q)-min(q)))
        elif new_x < min(q):
#            print 'Wrapping {} to {}'.format(new_x, max(q) -
#                                             abs((new_x-min(q))))
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
# HISTOGRAM POSITIONS
##################################


##################################
# COUNT FLUX
##################################

##################################
# PLOTTING
##################################
'''
fig, gs, axes = generate_axes_pad(nrows=num_landscapes, ncols=4, v_pad=0.3,
                                  h_pad=0.2, figsize=(12, 6))
for i in range(num_landscapes):
    if i % 2 == 0:
        c = 'r'
        t = 'Apo'
    else:
        c = 'b'
        t = 'Holo'
    ax = axes[i][0]
    ax.plot(q, E[i], color=c, lw=2)
    ax.set_title(t + ' energy landscape')
    ax = axes[i][1]
    lbl = u'$x(t+1) = x(t) + \\frac{D}{kT} F \\Delta t + \\Gamma$'
    ax.scatter(positions, zeros(len(positions)), label=lbl, color='k',
               alpha=0.01)
#    ax.plot(range(len(positions)), positions, color='k')
    ax.set_ylim(-1, 1)
    ax.set_title('Positions along reaction coordinate')
    ax.legend()
    ax = axes[i][2]
    counts, edges = histogram(positions, range=(min(q), max(q)),
                              bins=10, normed=True)
    mids = (edges[1:]+edges[:-1])/2.
    ax.bar(mids, counts, color=c, width=mids[1]-mids[0])
    ax.set_title('Histogram of positions')

    ax = axes[i][3]
    ax.plot(positions[::-1], range(len(positions))[::-1], color='k')
    ax.set_ylabel('Time')
    ax.set_xlabel('Position')
    ax.set_title('Kymograph')
ax.set_xlabel('Arbitrary reaction coordinate')
plt.show()
'''
c='r'


#fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
#                                  h_pad=0.2, figsize=(12, 12))
#ax = axes[0][0]
#ax.plot(range(len(boltzmann)), [i/sum(boltzmann) for i in boltzmann],
#        color='k')
#plt.show()

fig, gs, axes = generate_axes_pad(nrows=2, ncols=2, v_pad=0.3,
                                  h_pad=0.2, figsize=(12, 12))
ax = axes[0][0]
ax.plot(q, E[0], color=c, lw=2)
ax.set_title('Energy landscape')

ax = axes[0][1]
# lbl = u'$x(t+1) = x(t) + \\frac{D}{kT} F \\Delta t + \\Gamma$'
# ax.scatter(positions, zeros(len(positions)), label=lbl, color='k',
#            alpha=0.01)
# #    ax.plot(range(len(positions)), positions, color='k')
# ax.set_ylim(-1, 1)
# ax.set_title('Positions along reaction coordinate')
# ax.legend()
ax.plot(q[1:], force, color=c, lw=2)
ax.set_title('Force')

ax = axes[1][0]
counts, edges = histogram(positions, range=(min(q), max(q)),
                          bins=10, normed=True)
mids = (edges[1:]+edges[:-1])/2.
ax.bar(mids, counts, color=c, width=mids[1]-mids[0])

reduced_bins = 11  # This is num bins + 1
bins = linspace(0, len(boltzmann), num=reduced_bins)
reduced = [sum(boltzmann[int(ceil(bins[i])):int(floor(bins[i+1]))])
           for i in range(len(bins)-1)]
offset = abs(mids[1]-mids[0])/2
ax.plot(mids+offset, [i/(sum(reduced)*(edges[1]-edges[0])) for i in reduced],
        color='k', alpha=0.5, marker='o', ls='--', label='Boltzmann')
#ax.bar(mids, reduced, color='k', alpha=0.5,
#       width=mids[1]-mids[0])
ax.set_title('Histogram of positions')
ax.legend()

ax = axes[1][1]
counts, edges = histogram(jumps, range=(min(jumps), max(jumps)),
                          bins=10, normed=True)
mids = (edges[1:]+edges[:-1])/2.
ax.bar(mids, counts, color=c, width=mids[1]-mids[0])
ax.set_title('Histogram of step sizes')
plt.show()

fig, gs, axes = generate_axes_pad(nrows=1, ncols=1, v_pad=0.3,
                                  h_pad=0.2, figsize=(12, 12))
ax = axes[0][0]
cmap = cm.jet
c = linspace(min(positions), max(positions), len(positions))
cax = ax.scatter(positions[::-1], range(len(positions))[::-1], c=c, cmap=cmap,
                 s=40, lw=0, alpha=0.2,
                 label='Net boundary crossings = {} \n Total boundary'
                 ' crossings = {}'
                 .format(sum(fluxes), sum(absolute(fluxes))))
ax.set_title('Kymograph')
ax.set_xlabel('Position')
ax.set_ylabel('Time (number of steps)')
ax.legend()
cbar = fig.colorbar(cax, ticks=[min(positions), max(positions)])
cbar.ax.set_yticklabels(['End', 'Start'])
plt.show()