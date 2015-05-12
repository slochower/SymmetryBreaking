def initialize_2DBD():
    print('Creating 2D energy surface and initializing simulation.')
    x = q
    y = r
    z_one = energy(x, shift=0)
    z_two = energy(x, shift=pi)
    # Sloppy overloading of the name energy... 
    z = zeros((len(y), len(z_one)))

    for j in y:
        percent = j/transition_distance
        z[where(y==j)[0][0]] = percent*z_one + (1-percent)*z_two

    ####################################################################
    # Deprecated
    #forces = [[force_2DBD(i, dx, z[0][:], j, dy, z[:][0]) 
    #          for i in range(len(x)-1)] for j in range(len(y)-1)]
    ####################################################################

    # This is more ugly, but almost as fast, and much easier to 
    # read if I separate out the forces into different arrays.
    x_forces = array([[force(x_vals, dx, z[y_vals]) for x_vals in range(len(q)-1)] for y_vals in range(len(r)-1)])
    y_forces = array([[force(y_vals, dy, z[:,x_vals]) for y_vals in range(len(r)-1)] for x_vals in range(len(q)-1)]).T

    # Need to transpose for y in imshow, but not for x... (?)
    #    In [108]: plt.figure; plt.imshow(array(x_forces),origin='lower left'); plt.xlabel('x (elements, not distance)'); plt.ylabel('y (elements, not distance)'); plt.title('Forces in the $\hat{x}$ direction'); plt.colorbar(); plt.show()
    # In [106]: plt.figure; plt.imshow(array(y_forces).T,origin='lower left'); plt.xlabel('x (elements, not distance)'); plt.ylabel('y (elements, not distance)'); plt.title('Forces in the $\hat{y}$ direction'); plt.colorbar(); plt.show()

    boltzmann = [exp(-z[i]/kT) for i in range(len(z))]
    pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(len(boltzmann))]

    walker, net_flux = zeros((1, total_timesteps)), zeros((1, total_timesteps))
    start = [0.0, 0.0]
    steps_executed = 0

    return(z, x_forces, y_forces, boltzmann, pdf, walker, net_flux, steps_executed, start)


def column(matrix, i):
    return [row[i] for row in matrix]


def force_lookup(q, array, value):
    push = interp(value, q, array)
    return(push)

#@profile
def force_lookup_2D(force_interpolation, x, y):
    print('Interpolating force')
    return force_interpolation(x,y)

#@profile    
def simulate_2DBD(x, dx, y, dy, D, kT, dt, x_forces, y_forces, energies,steps_on_this_landscape, **kwargs):

    print('Calling simulation')
    # Set up force interpolation functions once...
    F_x = interpolate.interp2d(q[:-1], r[:-1], x_forces)
    F_y = interpolate.interp2d(q[:-1], r[:-1], y_forces)


    positions = []
    positions.append(hstack((x, y)))
    t = 1
    while t < steps_on_this_landscape-1:
        print(t)
        #######################################
        # BROWNIAN WALK
        #######################################
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        #F = force_lookup_2D(F_x, x, y)
        F = F_x(x,y)
        new_x = x + (D / kT) * F * dt + g
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        #F = force_lookup_2D(F_y, x, y)
        F = F_y(x,y)
        new_y = y + (D / kT) * F * dt + g

        # BOUNDARY CONDITIONS NEED TO BE HANDLED.

        ####################
        # RECORD KEEPING
        ####################
        t += 1
        positions.append(hstack((new_x, new_y)))
        x = new_x
        y = new_y

    return(positions)




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

