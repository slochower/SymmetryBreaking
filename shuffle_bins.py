#!/usr/bin/python

from math import *
from numpy import *
from scipy import optimize

execfile('../common/figures.py')


def landscapes(x, shift=0, noise=0):
    return float(cos(2*x + shift * pi) + random.normal(0, 0.05 * pi, 1))

#################################
# PARAMETERS
#################################
num_landscapes = 4
x = arange(0, 2 * pi, 0.01)
shift = pi                # Shift of landscape between holo and apo

##################################
# BUILD ENERGY LANDSCAPES
##################################
E = []    # Container for energy landscapes
p_x = []  # Container for optimization output parameters (x values)
p_y = []  # Container for optimization output parameters (y values)
s = []    # List of beginning points for moving on the energy landscapes

for i in range(num_landscapes):
    if i % 2 == 0:
        j = 0.0
    else:
        j = shift
    if i == 0:
        start = 0.1
    else:
        start = p_x[i - 1]
    E.append([landscapes(i, shift=j) for i in x])
    p = optimize.minimize(landscapes, start, method='Powell',
                          args=(j,))
    '''
    See: http://en.wikipedia.org/wiki/Powell%27s_method
    The method minimises the function by a bi-directional search along each
     search vector, in turn.

    The method is useful for calculating the local minimum of a continuous but
     complex function, especially one without an underlying mathematical
     definition, because it is not necessary to take derivatives.
    '''
    p_x.append(p['x'])
    p_y.append(p['fun'])
    s.append(start)

##################################
# NORMING AND BINNING
##################################
F = [[j/sum(E[i]) for j in E[i]] for i in range(num_landscapes)]

#    counts, edges = histogram(E[i], range=(0, max(x)), bins=100, normed=True)
#    mids = (edges[:-1]+edges[1:])/2.

##################################
# PLOTTING
##################################
fig, gs, axes = generate_axes_pad(nrows=num_landscapes, ncols=2, v_pad=0.3,
                                  h_pad=0.2, figsize=(6, 24))
for i in range(num_landscapes):
    if i % 2 == 0:
        c = 'r'
        t = 'Apo'
    else:
        c = 'b'
        t = 'Holo'
    ax = axes[i][0]
    ax.plot(x, E[i], color=c, lw=2)
    ax.scatter(p_x[i], p_y[i], color='k', s=40, zorder=10)
    ax.axvline(x=s[i], color='k', ls='dashed')
    ax.text(s[i]+0.1, 1.5, ' Start')
    ax.set_title(t + ' energy landscape')
    ax = axes[i][1]
    ax.plot(x, F[i], color=c, lw=2)
ax.set_xlabel('Arbitrary reaction coordinate')
plt.show()
