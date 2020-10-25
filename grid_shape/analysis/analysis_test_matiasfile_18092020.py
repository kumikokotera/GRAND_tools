import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

path = "/Users/benoitl/Documents/GRAND/Data_grids/20200918/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/Trihex"


threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 
 
### "Creating" the layout
pos, offset = grids.create_grid_univ("trihex", 125, do_prune=False, input_n_ring=10)

## Creating the pruned layout
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 125, do_prune=True, input_n_ring=10)

# creating the all mask (all the antennas)
mask  = masks.make_all_mask(input_n_ring=10)


lay1 = layout.Layout(path, pos, mask, "all", threshold, n_trig_thres)
lay2 = layout.Layout(path, pos2, mask2, "simple", threshold, n_trig_thres)

# creating a random mask with only 5% of the n_ring=10 antennas kept
mask_rand_5 = masks.make_mask_random(input_n_ring=10, n_keep_ratio=0.05)

# creating the layout associated with this mask
lay_rand_5 = layout.Layout(path, pos, mask_rand_5, "rand_5", threshold, n_trig_thres)

# creating the trihex 250 grid out of the trixhex 125
mask_tri250 = masks.make_trihex_new_out_of_125(pos, 250, 10) 
lay_tri250 = layout.Layout(path, pos, mask_tri250, "tri250", threshold, n_trig_thres)

lay_tri250.plot_layout()

# creating the trihex 500 grid out of the trixhex 125
mask_tri500 = masks.make_trihex_new_out_of_125(pos, 500, 10) 
lay_tri500 = layout.Layout(path, pos, mask_tri500, "tri500", threshold, n_trig_thres)

lay_tri500.plot_layout()



# creating the trihex 500 grid out of the trixhex 125
mask_tri1000 = masks.make_trihex_new_out_of_125(pos, 1000, 10) 
lay_tri1000 = layout.Layout(path, pos, mask_tri1000, "tri1000", threshold, n_trig_thres)

lay_tri1000.plot_layout()







plt.figure(3, figsize=(8,6))
plt.clf()
n_zenith_bins = len(lay1.zenith_bins_centers)
n_bins_to_plot = 4
i_bins = [(i) * n_zenith_bins // (n_bins_to_plot+1) for i in range(n_bins_to_plot)]
#i_bins = [0, 1, 2, 3, 4]
for k, i_bin in enumerate(i_bins):
    plt.errorbar(
        np.log10(lay1.energy_bins_centers*1e18),
        lay1.detection_rate[:,i_bin],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=8,
        color = "C1",
        ls = '-',
        label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[i_bin],lay1.zenith_bins_limits[i_bin+1])
    )
    plt.errorbar(
        np.log10(lay1.energy_bins_centers*1e18),
        lay_rand_5.detection_rate[:,i_bin],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=8,
        color="C2",
        ls = '-',
        # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
    )
    plt.errorbar(
        np.log10(lay1.energy_bins_centers*1e18),
        lay_tri250.detection_rate[:,i_bin],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=8,
        color="C3",
        ls = '-',
        # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
    )
    plt.errorbar(
        np.log10(lay1.energy_bins_centers*1e18),
        lay_tri500.detection_rate[:,i_bin],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=8,
        color="C4",
        ls = '-',
        # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
    )
    plt.errorbar(
        np.log10(lay1.energy_bins_centers*1e18),
        lay_tri1000.detection_rate[:,i_bin],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=8,
        color="C5",
        ls = '-',
        # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
    )

    plt.yscale('log')
    plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} km^{-2}]$')
    plt.xlabel('log E/eV')
    legend1 = plt.legend(loc=1)
    #splt.xscale('log')

custom_lines = [
    Line2D([0], [0], color="C1", lw=4),
    Line2D([0], [0], color="C2", lw=4),
    Line2D([0], [0], color="C3", lw=4),
    Line2D([0], [0], color="C4", lw=4),
    Line2D([0], [0], color="C5", lw=4)
]

legend2 = plt.legend(custom_lines, ['trihex all', 'random 5%',"trihex 250", "trihex500", "trihex1000"], loc=3)

plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)





plt.figure(4, figsize=(8,6))
plt.clf()
n_energy_bins = len(lay1.energy_bins_centers)
n_bins_to_plot = 4
i_bins = [(i+1) * n_energy_bins // (n_bins_to_plot+1) for i in range(n_bins_to_plot)]
#i_bins = [0, 1, 2, 3, 4]
for k, i_bin in enumerate(i_bins):
#for iener in range(0, len(lay1.energy_bins_limits)-1):
    plt.errorbar(
        lay1.zenith_bins_centers,
        lay1.detection_rate[i_bin,:],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=10,
        ls = '-',
        color="C1",
        label='%4.2f > log E/eV >%4.2f'%(
            np.log10(lay1.energy_bins_limits[i_bin]*1e18),
            np.log10(lay1.energy_bins_limits[i_bin+1]*1e18)
        )
    )
    plt.errorbar(
        lay_rand_5.zenith_bins_centers,
        lay_rand_5.detection_rate[i_bin,:],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=10,
        ls = '-',
        color="C2"
        )
    plt.errorbar(
        lay_tri250.zenith_bins_centers,
        lay_tri250.detection_rate[i_bin,:],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=10,
        ls = '-',
        color="C3"
    )
    plt.errorbar(
        lay_tri500.zenith_bins_centers,
        lay_tri500.detection_rate[i_bin,:],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=10,
        ls = '-',
        color="C4"
    )
    plt.errorbar(
        lay_tri1000.zenith_bins_centers,
        lay_tri1000.detection_rate[i_bin,:],
        fmt=SYM_LIST[i_bin],
        capsize=2,
        alpha=0.7,
        ms=10,
        ls = '-',
        color="C5"
    )


    plt.yscale('log')
    plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} km^{-2}]$')
    plt.xlabel('Zenith [deg]')
    #plt.title('%s, %4.2f > zenith > %4.2f deg'%(lay1.mask_name, lay1.[izen], thetal[izen+1]))
    legend1 = plt.legend(loc=1)
    
 
custom_lines = [
    Line2D([0], [0], color="C1", lw=4),
    Line2D([0], [0], color="C2", lw=4),
    Line2D([0], [0], color="C3", lw=4),
    Line2D([0], [0], color="C4", lw=4),
    Line2D([0], [0], color="C5", lw=4)
]

legend2 = plt.legend(custom_lines, ['trihex all', 'random 5%',"trihex 250", "trihex500", "trihex1000"], loc=3)

plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)



 ## make a trixhex 250 nring =10 lauout


