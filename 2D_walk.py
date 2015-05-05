from mpl_toolkits.mplot3d.axes3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X = arange(0, 2*pi, 0.25)
Y = arange(0, 40, 1)
x, y = meshgrid(X, Y)
z = landscapes(x)
z[-1] = landscapes(x[-1], shift=pi)
z[1:-1] = [mean(energy[0]) for k in range(shape(z)[1])]
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
from scipy.interpolate import griddata
grid_x, grid_y = mgrid[0:2*pi:0.25, 0:41:1]

zero = array([q, [0]*629, landscapes(q)]).T
one = array([q, [40]*629, landscapes(q, shift=pi)]).T
values = vstack((zero, one))
points = tuple([tuple(row) for row in values[:,0:2]])
grid_z0 = griddata(array(points), values[:,2], (grid_x, grid_y), method='nearest')

surf = ax.plot_surface(grid_x, grid_y, grid_z0, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)
plt.show()    


# 1. Create a 3D mesh.
# 2. Assign some value to some points in the mesh.
# 3. Interpolate (?)

# or

# 1. Define two curves.
# 2. Find way to morph between them.
# 3. Space out the two curves.