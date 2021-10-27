import matplotlib.pyplot as plt
import os

from grid_shape_lib.modules import grids as grids
from grid_shape_lib.modules import masks as masks


def plot_mask(mask, nring, plot_suffix):
    pos, _ = grids.create_grid_univ("trihex", 250, input_n_ring=nring)
    plt.figure()
    plt.scatter(
        pos[0],
        pos[1],
        c = 1-mask[:,0],
        s = 3
    )
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.savefig(
        os.path.join("layout_{}.png".format(plot_suffix))
    )


nring = 10
# All mask
mask1 = masks.make_all_mask(nring)
plot_mask(mask1, 10, 'all_mask')

# mask_random
mask2 = masks.make_mask_random(10, .1)
plot_mask(mask2, 10, 'random10')

#make_coarse_out_of_fine_trihex
mask3 = masks.make_coarse_out_of_fine_trihex(250, 'hexhex', 125, 20)
plot_mask(mask3, 20, 'make_coarse_out_of_fine_trihex')


#make_coarse_out_of_fine_trihex2
mask3 = masks.make_coarse_out_of_fine_trihex_2(250, 'hexhex',5 , 125, 20)
plot_mask(mask3, 20, 'make_coarse_out_of_fine_trihex2')


# make_mini_islands
# (islands are located at positions of a trihexstep1000grid)
mask4 = masks.make_mini_islands(
    125, 20, with_center=False)
plot_mask(mask4, 20, 'make_make_mini_islands')

#make_limited_mini_island
mask5 = masks.make_limited_mini_island(
    'hexhex',
    250,
    1000,
    "trihex",
    "trihex",
    5,
    250,
    125, 
    20
)
plot_mask(mask5, 20, 'make_limited_mini_islands')


#full trihex nring =0, 1,2, step 125 out of a finer nring=10 125 trihex

mask6 = masks.make_coarse_out_of_fine_trihex_2(125, 'trihex', 1, 125, 10)
plot_mask(mask6, 10, 'make6')