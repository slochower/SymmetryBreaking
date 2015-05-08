def initialize_2DBD():
    x = q
    y = r
    z_one = energy(x, shift=0)
    z_two = energy(x, shift=pi)
    # Sloppy overloading of the name energy... 
    z = zeros((len(y), len(z_one)))

    for j in y:
        percent = j/transition_distance
        z[where(y==j)[0][0]] = percent*z_one + (1-percent)*z_two

    force = [[force_2DBD(i, dx, z[0][:], j, dy, z[:][0]) 
              for i in range(len(x)-1)] for j in range(len(y)-1)]
    boltzmann = [exp(-z[i]/kT) for i in range(len(z))]
    pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(len(boltzmann))]

    steps_executed, landscapes_sampled, start = 0, 0, 0.0
    steps_on_A, steps_on_B = 0, 0
    # Hoping I can store tuples in walker[] and have it be okay...
    walker, net_flux = zeros((1, total_timesteps)), zeros((1, total_timesteps))

    start = (0.0, 0.0)

    return(z, force, boltzmann, pdf, walker, net_flux, steps_executed, landscapes_sampled, start, steps_on_A, steps_on_B)


def force_2DBD(i, dx, x_array, j, dy, y_array):
    F_x = -((x_array[i+1] - x_array[i]) / dx)
    F_y = -((y_array[j+1] - y_array[j]) / dy)
    return(F_x, F_y)


def force_lookup(q, array, value):
    push = interp(value, q, array)
    return(push)

    
def simulate_2DBD(x, dx, y, dy, D, kT, dt, forces, energies,steps_on_this_landscape, **kwargs):

    positions = []
    positions.append((x, y))
    t = 1
    while t < steps_on_this_landscape-1:
        #######################################
        # BROWNIAN WALK
        #######################################
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        # x forces at y = 0, I think... 
        # No, I don't think is correct anymore.
        # This is a roadblock for the simulation.
        # In [48]: shape(forces)
        # Out[48]: (99, 628, 2)


        x_forces = [x[0] for x in force[0]]
        F = force_lookup(q[:-1], force, x)
        new_x = x + (D / kT) * F * dt + g
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(r[:-1], force, y)
        new_y = y + (D / kT) * F * dt + g
        ####################
        # RECORD KEEPING
        ####################
        t += 1
        fluxes.append(swap)
        positions.append(new_x, new_y)
        x = new_x
        y = new_y

    return(positions, fluxes)




if 'plotting' in options:
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt 

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = meshgrid(x_transition, y_transition)
    surf = ax.plot_surface(X, Y, z_transition, rstride=1, cstride=1, cmap=cm.coolwarm,  linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    from mayavi import mlab
    fig = mlab.figure(size=(2000,2000))
    plot = mlab.mesh(X, Y, z_transition)
    # mlab.mesh(X, Y, z_transition, representation='wireframe', color=(0, 0, 0))
    mlab.axes(plot)
    mlab.colorbar(plot, title='Energy', orientation='vertical')

    duration = 10 
    def make_frame(t):
        """ Generates and returns the frame for time t. """
        mlab.view(azimuth=(360*t/duration), distance=85) # camera angle
        return mlab.screenshot(antialiased=True) # return a RGB image

    import  moviepy.editor as mpy
    animation = mpy.VideoClip(make_frame, duration=duration).resize(0.5)
    animation.write_videofile("2D_rotation_grey.mp4", fps=20)

    plot = mlab.mesh(X, Y, z_transition)
    # mlab.mesh(X, Y, z_transition, representation='wireframe', color=(0, 0, 0))
    mlab.axes(plot)
    mlab.colorbar(plot, title='Energy', orientation='vertical')

