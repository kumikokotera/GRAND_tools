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


### Commented analysis script

path = "/Users/benoitl/Documents/GRAND/Data_grids/3_dec_20/"
plot_path = '/Users/benoitl/Documents/GRAND/Data_grids/3_dec_20/plots'
os.makedirs(plot_path, exist_ok=True)
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/Trihex"

## use 75 !!!!
threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 
 
### "Creating" the layout
pos, offset = grids.create_grid_univ("trihex", 250, do_prune=False, input_n_ring=20)

## Creating the pruned layout
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 250, do_prune=True, input_n_ring=20)


# creating the all mask (all the antennas)
mask  = masks.make_all_mask(input_n_ring=20)

mask_mini_islands = masks.make_mini_islands(250, 20)


mask_1ring_250 = masks.make_coarse_out_of_fine_trihex_2(250, 'hexhex', 0, 250, 20)
mask_1ring_500 = masks.make_coarse_out_of_fine_trihex_2(500, 'hexhex', 0, 250, 20)
mask_1ring_1000 = masks.make_coarse_out_of_fine_trihex_2(1000, 'hexhex', 0, 250, 20)



mask_limited_mini_island = masks.make_limited_mini_island(
    'hexhex',
    250,
    1000,'hexhex',"trihex", 6, 250, 250, 20)

mask_tri1000 = masks.make_coarse_out_of_fine_trihex(1000, "trihex", 250, 20)
mask_hex1000 = masks.make_coarse_out_of_fine_trihex(1000, "hexhex", 250, 20)
mask_hex500 = masks.make_coarse_out_of_fine_trihex(500, "hexhex", 250, 20)
mask_hex250 = masks.make_coarse_out_of_fine_trihex(250, "hexhex", 250, 20)


mask_hex1250 = masks.make_coarse_out_of_fine_trihex(1250, "hexhex", 250, 20)


mask_hex1000_island250 =  masks.make_central_island_out_of_fine(
    1000,
    "hexhex",
    250,
    "hexhex",
    3,
    250,
    20
)


mask_hex1000_triisland250 =  masks.make_central_island_out_of_fine(
    1000,
    "hexhex",
    250,
    "trihex",
    2,
    250,
    20
)


mask_hex1250_triisland250 =  masks.make_central_island_out_of_fine(
    1250,
    "hexhex",
    250,
    "trihex",
    3,
    250,
    20
)


lay_all = layout.Layout(path, pos, mask, "all", threshold, n_trig_thres, input_n_ring = 20, primary = "Proton")
    


lay_limited_mini_island = layout.Layout(
    path,
    pos,
    mask_limited_mini_island,
    "lay_limited_mini_island",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton"
)
lay1 = lay_limited_mini_island

if False:
    lay_hex1000_island250 = layout.Layout(
        path,
        pos,
        mask_hex1000_island250,
        "lay_hex1000_island250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay.plot_2d_efficiency(7, 7)

if False:
    print('doing: lay_1ring_250')
    lay_1ring_250 = layout.Layout(
        path,
        pos,
        mask_1ring_250,
        "lay_1ring_250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_1ring_250.make_all_plots()
    lay_1ring_250.plot_2d_efficiency(7,7)

if False:
    print('doing: lay_1ring_500')
    lay_1ring_500 = layout.Layout(
        path,
        pos,
        mask_1ring_500,
        "lay_1ring_500",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_1ring_500.make_all_plots()
    lay_1ring_500.plot_2d_efficiency(7, 7)

if False:
    print('doing: lay_1ring_1000')
    lay_1ring_1000 = layout.Layout(
        path,
        pos,
        mask_1ring_1000,
        "lay_1ring_1000",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_1ring_1000.make_all_plots()
    lay_1ring_1000.plot_2d_efficiency(7, 7)


if False:
    print('doing: lay_hex1000_triisland250')
    lay_hex1000_triisland250 = layout.Layout(
        path,
        pos,
        mask_hex1000_triisland250,
        "lay_hex1000_triisland250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )

    lay_hex1000_triisland250.make_all_plots()
    lay_hex1000_triisland250.plot_2d_efficiency(7, 7)



if True:
    print('doing: lay_hex1250_triisland250')
    lay_hex1250_triisland250 = layout.Layout(
        path,
        pos,
        mask_hex1250_triisland250,
        "lay_hex1250_triisland250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_hex1250_triisland250.make_all_plots()
    lay_hex1250_triisland250.plot_2d_efficiency(7, 7)

if False:
    lay_tri1000 = layout.Layout(
        path,
        pos,
        mask_tri1000,
        "lay_tri1000",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )

if False:
    lay_hex1000 = layout.Layout(
        path,
        pos,
        mask_hex1000,
        "lay_hex1000",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_hex1000.make_all_plots()
    lay_hex1000.plot_2d_efficiency(7, 7)


    lay_hex500 = layout.Layout(
        path,
        pos,
        mask_hex500,
        "lay_hex500",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_hex500.make_all_plots()
    lay_hex500.plot_2d_efficiency(7, 7)


if False:
    lay_hex1250 = layout.Layout(
        path,
        pos,
        mask_hex1250,
        "lay_hex1250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_hex1250.make_all_plots()

    lay1 = layout.Layout(path, pos, mask, "all", threshold, n_trig_thres, input_n_ring = 20, primary = "Proton")
    lay2 = layout.Layout(path, pos2, mask2, "simple", threshold, n_trig_thres,  input_n_ring = 20, primary = "Proton")

    lay1.plot_2D_differential_rate()



    lay_mini_islands2 = layout.Layout(
        path,
        pos,
        mask_mini_islands,
        "mini_islands2",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )



    mask_spiral = masks.make_spiral_mask(1000, 1.4)
    lay_spiral = layout.Layout(
        path,
        pos,
        mask_spiral,
        "spiral",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_spiral.make_all_plots()



    mask_spiral2 = masks.make_spiral_mask(1000, 1.05)
    lay_spiral2 = layout.Layout(
        path,
        pos,
        mask_spiral2,
        "spiral2",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_spiral2.make_all_plots()



    mask_spiral3 = masks.make_spiral_mask_with_minipose(1000, 1.7, 21)
    lay_spiral3 = layout.Layout(
        path,
        pos,
        mask_spiral3,
        "spiral3",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_spiral3.make_all_plots()

if True:
    mask_spiral4 = masks.make_spiral_mask_with_minipose(1000, 1.8, 21)
    lay_spiral4 = layout.Layout(
        path,
        pos,
        mask_spiral4,
        "spiral4",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay_spiral4.make_all_plots()

if False:
    mask_spiral5 = masks.make_spiral_mask_with_minipose(1000,2.2, 17, "trihex")
    lay_spiral5 = layout.Layout(
        path,
        pos,
        mask_spiral5,
        "spiral5",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )


    lay_spiral5.make_all_plots()


    lay_mini_islands2.make_all_plots()

    lay_tri1000.make_all_plots()


    lay_hex1000.make_all_plots()

    lay1.make_all_plots()

    lay_hex1000_island250.make_all_plots()


    lay_limited_mini_island.make_all_plots()

