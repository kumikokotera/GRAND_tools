import hexy as hx
import numpy as np
import matplotlib.pyplot as plt
import grids as grids
radius = 1000  # [m]

radius = 1291.472  # [m]

area = hx.get_area(radius) * 3
cube0 = np.array(
    [
        [0, 0, 0],
        [1, 0, -1],
        [0, 1, -1],
    ]
)

pos = hx.cube_to_pixel(cube0, radius)
corners = hx.get_corners(pos, radius)


sh = np.array(corners).shape
corners = corners.transpose(0, 2, 1)
corners = np.array(corners).reshape((sh[0]*sh[2], sh[1]))

corners_x, corners_y = grids.remove_redundant_point(corners[:, 0], corners[:, 1])

pos_gp13 = np.array([corners_x, corners_y]).transpose()

plt.figure(1)
plt.clf()
plt.plot(pos[:, 0], pos[:, 1], 'k.', label='hex. center')
plt.plot(corners_x, corners_y, 'r.', label='antenna position')
plt.axis('equal')
plt.legend(loc=0)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()
plt.title('radius = {} m '.format(radius))
plt.savefig('layout_gp13.png')


np.savetxt('pos_layout_gp13.txt', pos_gp13)

