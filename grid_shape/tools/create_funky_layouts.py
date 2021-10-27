import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
from grid_shape import layout as layout
from grid_shape import masks as masks





SYM_LIST = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
MYC = ['0','0.20','0.4','0.6','0.8']



pos, offset = grids.create_grid_univ("trihex", 125, do_prune=False, input_n_ring=10)

## Creating the pruned layout
##pos2, offset2, mask2 = grids.create_grid_univ("trihex", 125, do_prune=True, input_n_ring=10)

# creating the all mask (all the antennas)
mask  = masks.make_all_mask(input_n_ring=10)



# creating the trihex 500 grid out of the trixhex 125
mask_tri1000 = masks.make_trihex_new_out_of_125(pos, 1000, 10) 


mask_island1 = masks.make_centralisland_out_of_125(pos, 10) 
mask_island2 = masks.make_centralisland_out_of_125_v2(pos, 10) 





minipos, _ = grids.create_grid_univ("trihex", step, do_prune=False, input_n_ring=0)
minipos2 = minipos.copy()
id1000 = np.where(mask_tri1000)[0]

for id in id1000:
    minipos2 = np.hstack([minipos2,  minipos + pos[:,id].reshape(3,1)])

#minipos2 = minipos + pos[:,id1000[0]].reshape(3,1)
#minipos3 = minipos + pos[:,id1000[1]].reshape(3,1)
#minipos4 = minipos + pos[:,id1000[2]].reshape(3,1)
#minipos5 = minipos + pos[:,id1000[3]].reshape(3,1)
#minipos6 = minipos + pos[:,id1000[4]].reshape(3,1)

#notsominipos = np.hstack([minipos2, minipos3, minipos4, minipos5, minipos6])

notsomini_mask = masks.prune_pos(pos, minipos2)


plt.figure(1,figsize=(8, 6))
plt.clf()
plt.scatter(
    pos[0],
    pos[1],
    c = 1-notsomini_mask[:,0],
    s = 3
)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')




plt.figure(2,figsize=(8, 6))
plt.clf()
plt.scatter(
    pos[0],
    pos[1],
    c = 1-mask_island1[:,0],
    s = 3
)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')