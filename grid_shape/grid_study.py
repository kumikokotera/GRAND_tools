import numpy as np
import matplotlib.pyplot as plt
import grids

n = 1000

steprect = 275
pos0, offset0 = grids.create_grid_univ('rect', steprect, angle=0, do_offset=False)

d0_list = []
d1_list = []
d2_list = []
d3_list = []
d4_list = []

for i in range(0,n):
    pos, offset = grids.create_grid_univ('rect', steprect, angle=0, do_offset=True)
    d = np.sqrt((pos0[0,:]-offset[0])**2 + (pos0[1,:]-offset[1])**2)

    ind = np.argsort(d)
    d0_list.append(d[ind[0]])
    d1_list.append(d[ind[1]])
    d2_list.append(d[ind[2]])
    d3_list.append(d[ind[3]])
    d4_list.append(d[ind[4]])

d0 = np.array(d0_list)
d1 = np.array(d1_list)
d2 = np.array(d2_list)
d3 = np.array(d3_list)
d4 = np.array(d4_list)


stephex = 250
pos0, offset0 = grids.create_grid_univ('hexhex', stephex, angle=0, do_offset=False)

d0hex_list = []
d1hex_list = []
d2hex_list = []
d3hex_list = []
d4hex_list = []

for i in range(0,n):
    pos, offset = grids.create_grid_univ('hexhex', stephex, angle=0, do_offset=True)
    d = np.sqrt((pos0[0,:]-offset[0])**2 + (pos0[1,:]-offset[1])**2)

    ind = np.argsort(d)
    d0hex_list.append(d[ind[0]])
    d1hex_list.append(d[ind[1]])
    d2hex_list.append(d[ind[2]])
    d3hex_list.append(d[ind[3]])
    d4hex_list.append(d[ind[4]])


d0hex = np.array(d0hex_list)
d1hex = np.array(d1hex_list)
d2hex = np.array(d2hex_list)
d3hex = np.array(d3hex_list)
d4hex = np.array(d4hex_list)




plt.figure(1) 
plt.clf()
plt.hist(d0, bins=20, range=[0,275], label='rect, step = %d m'%(np.int32(steprect)))
plt.hist(d0hex, bins=20, range=[0,275], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
plt.xlabel('antenna distance to core [m]')
plt.ylabel('N')
plt.legend()
plt.show()

plt.figure(2) 
plt.clf()
plt.hist(d1, bins=20, range=[0,275], label='rect, step = %d m'%(np.int32(steprect)))
plt.hist(d1hex, bins=20, range=[0,275], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
plt.xlabel('antenna distance to core [m]')
plt.ylabel('N')
plt.legend()
plt.show()

plt.figure(3) 
plt.clf()
#plt.hist(d2, bins=20, range=[0,500], label='rect, step = %d m'%(np.int32(steprect)))
#plt.hist(d2hex, bins=20, range=[0,500], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
plt.hist(d2, bins=20, range=[0,500], label='rect, step = %d m'%(np.int32(steprect)))
plt.hist(d2hex, bins=20, range=[0,500], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
plt.xlabel('antenna distance to core [m]')
plt.ylabel('N')
plt.legend()
plt.show()

plt.figure(4) 
plt.clf()
#plt.hist(d2, bins=20, range=[0,500], label='rect, step = %d m'%(np.int32(steprect)))
#plt.hist(d2hex, bins=20, range=[0,500], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
#plt.hist(d3, bins=20, range=[0,500], label='rect, step = %d m'%(np.int32(steprect)))
#plt.hist(d3hex, bins=20, range=[0,500], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
plt.hist(d4, bins=20, range=[0,500], label='rect, step = %d m'%(np.int32(steprect)))
plt.hist(d4hex, bins=20, range=[0,500], alpha=0.5, label='hex, step = %d m'%(np.int32(stephex)))
plt.xlabel('antenna distance to core [m]')
plt.ylabel('N')
plt.legend()
plt.show()