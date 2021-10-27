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


threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 
 
### "Creating" the layout
pos, offset = grids.create_grid_univ("trihex", 250, do_prune=False, input_n_ring=20)

## Creating the pruned layout
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 250, do_prune=True, input_n_ring=20)


# creating the all mask (all the antennas)
mask  = masks.make_all_mask(input_n_ring=20)

mask_mini_islands = masks.make_mini_islands(250, 20)

mask_limited_mini_island = masks.make_limited_mini_island(
    'hexhex',
    250,
    1000,'hexhex',"trihex", 6, 250, 250, 20)

mask_tri1000 = masks.make_coarse_out_of_fine_trihex(1000, "trihex", 250, 20)
mask_hex1000 = masks.make_coarse_out_of_fine_trihex(1000, "hexhex", 250, 20)
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

mask_hex3000_island250 =  masks.make_central_island_out_of_fine(
    3000,
    "hexhex",
    250,
    "hexhex",
    3,
    250,
    20
)

mask_hex6000_island250 =  masks.make_central_island_out_of_fine(
    6000,
    "hexhex",
    250,
    "hexhex",
    3,
    250,
    20
)

mask_hex6000_island250_2 =  masks.make_central_island_out_of_fine(
    6000,
    "hexhex",
    250,
    "hexhex",
    1,
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
if True:
    lay_hex3000_island250 = layout.Layout(
        path,
        pos,
        mask_hex3000_island250,
        "lay_hex3000_island250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay3 = lay_hex3000_island250
    lay_hex3000_island250.make_all_plots()

    lay_hex6000_island250 = layout.Layout(
        path,
        pos,
        mask_hex6000_island250,
        "lay_hex6000_island250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay6000 = lay_hex6000_island250
    lay_hex6000_island250.make_all_plots()

    lay_hex6000_island250_2 = layout.Layout(
        path,
        pos,
        mask_hex6000_island250_2,
        "lay_hex6000_island250_2",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton"
    )
    lay6000_2 = lay_hex6000_island250_2
    lay_hex6000_island250_2.make_all_plots()

if True:
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
    lay2 = lay_hex1000_triisland250
if False:

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




    # # creating a random mask with only 5% of the n_ring=10 antennas kept
    # mask_rand_5 = masks.make_mask_random(input_n_ring=20, n_keep_ratio=0.05)

    # # creating the layout associated with this mask
    # lay_rand_5 = layout.Layout(path, pos, mask_rand_5, "rand_5", threshold, n_trig_thres)

    # # creating the trihex 250 grid out of the trixhex 125
    # mask_tri250 = masks.make_trihex_new_out_of_125(pos, 250, 20) 
    # lay_tri250 = layout.Layout(path, pos, mask_tri250, "tri250", threshold, n_trig_thres)

    # lay_tri250.plot_layout()

    # # creating the trihex 500 grid out of the trixhex 125
    # mask_tri500 = masks.make_trihex_new_out_of_125(pos, 500, 20) 
    # lay_tri500 = layout.Layout(path, pos, mask_tri500, "tri500", threshold, n_trig_thres)

    # lay_tri500.plot_layout()



    # # creating the trihex 500 grid out of the trixhex 125
    # mask_tri1000 = masks.make_trihex_new_out_of_125(pos, 1000, 20) 
    # lay_tri1000 = layout.Layout(path, pos, mask_tri1000, "tri1000", threshold, n_trig_thres)

    # lay_tri1000.plot_layout()


    # mask_island1 = masks.make_centralisland_out_of_125(pos, 20) 
    # lay_island1 = layout.Layout(path, pos, mask_island1, "island1", threshold, n_trig_thres)

    # lay_island1.plot_layout()



    # mask_island2 = masks.make_centralisland_out_of_125_v2(pos, 10) 
    # lay_island2 = layout.Layout(path, pos, mask_island2, "island2", threshold, n_trig_thres)

    # lay_island2.plot_layout()








    plt.figure(3, figsize=(8,6))
    plt.clf()
    n_zenith_bins = len(lay1.zenith_bins_centers)
    n_bins_to_plot = 6
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
            lay2.detection_rate[:,i_bin],
            fmt=SYM_LIST[i_bin],
            capsize=2,
            alpha=0.7,
            ms=8,
            color = "C2",
            ls = '-',
            #label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[i_bin],lay1.zenith_bins_limits[i_bin+1])
        )
        # plt.errorbar(
        #     np.log10(lay1.energy_bins_centers*1e18),
        #     lay_rand_5.detection_rate[:,i_bin],
        #     fmt=SYM_LIST[i_bin],
        #     capsize=2,
        #     alpha=0.7,
        #     ms=8,
        #     color="C2",
        #     ls = '-',
        #     # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
        # )
        # plt.errorbar(
        #     np.log10(lay1.energy_bins_centers*1e18),
        #     lay_tri250.detection_rate[:,i_bin],
        #     fmt=SYM_LIST[i_bin],
        #     capsize=2,
        #     alpha=0.7,
        #     ms=8,
        #     color="C3",
        #     ls = '-',
        #     # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
        # )
        # plt.errorbar(
        #     np.log10(lay1.energy_bins_centers*1e18),
        #     lay_tri500.detection_rate[:,i_bin],
        #     fmt=SYM_LIST[i_bin],
        #     capsize=2,
        #     alpha=0.7,
        #     ms=8,
        #     color="C4",
        #     ls = '-',
        #     # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
        # )
        # plt.errorbar(
        #     np.log10(lay1.energy_bins_centers*1e18),
        #     lay_tri1000.detection_rate[:,i_bin],
        #     fmt=SYM_LIST[i_bin],
        #     capsize=2,
        #     alpha=0.7,
        #     ms=8,
        #     color="C5",
        #     ls = '-',
        #     # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
        # )

        plt.yscale('log')
        plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} km^{-2}]$')
        plt.xlabel('log E/eV')
        legend1 = plt.legend(loc=1)
        #splt.xscale('log')

    custom_lines = [
        Line2D([0], [0], color="C1", lw=4),
        Line2D([0], [0], color="C2", lw=4),
        #Line2D([0], [0], color="C3", lw=4),
        #Line2D([0], [0], color="C4", lw=4),
        #Line2D([0], [0], color="C5", lw=4)
    ]

    legend2 = plt.legend(
        custom_lines,
        ['trihex all', 'hexhex all'],#,"trihex 250", "trihex500", "trihex1000"],
        loc=3)

    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    plt.savefig(os.path.join(lay1.plot_path, "detection_rate_vs_energy.png"))




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
    plt.savefig(os.path.join(lay1.plot_path, "detection_rate_vs_zenith.png"))


    plt.figure(45)
    lay1.plot_layout()
    plt.savefig(os.path.join(lay1.plot_path, "lay1.png"))



    plt.figure(46)
    lay_rand_5.plot_layout(fig=46)
    plt.savefig(os.path.join(lay1.plot_path, "lay_rand_5.png"))


    plt.figure(47)
    lay_tri250.plot_layout(fig=47)
    plt.savefig(os.path.join(lay1.plot_path, "lay_tri250.png"))


    plt.figure(48)
    lay_tri500.plot_layout(fig=48)
    plt.savefig(os.path.join(lay1.plot_path, "lay_tri500.png"))



    plt.figure(49)
    lay_tri1000.plot_layout(fig=49)
    plt.savefig(os.path.join(lay1.plot_path, "lay_tri1000.png"))


    plt.figure(50)
    lay_island1.plot_layout(fig=50)
    plt.savefig(os.path.join(lay1.plot_path, "lay_island1.png"))


    plt.figure(51)
    lay_island2.plot_layout(fig=51)
    plt.savefig(os.path.join(lay1.plot_path, "lay_island2.png"))



    lay1.make_diff_event_rate(200)
    lay_tri1000.make_diff_event_rate(200)
    lay_tri500.make_diff_event_rate(200)
    lay_tri250.make_diff_event_rate(200)


    plt.figure(345)
    plt.clf()
    plt.plot(lay1.energy_bins_centers, lay1.integrated_ev_rate_no_trig, label="trig efficiency = 1")
    plt.plot(lay1.energy_bins_centers, lay1.integrated_ev_rate, label="trihex all")
    plt.plot(lay_tri250.energy_bins_centers, lay_tri250.integrated_ev_rate, label="trihex250")
    plt.plot(lay_tri500.energy_bins_centers, lay_tri500.integrated_ev_rate, label="trihex500")
    plt.plot(lay_tri1000.energy_bins_centers, lay_tri1000.integrated_ev_rate, label="trihex1000")
    plt.ylabel('Differential event rate [day'+"$^{-1}$"+"PeV"+"$^{-1}$"+']')
    plt.title('Detector area 200 km'+"$^2$")
    plt.xlabel('Energy [EeV]')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.savefig('diff_rate_200km2.png')




    ####   hist 2d

    plt.figure(100)
    plt.clf()
    plt.pcolor(
        np.log10(1e18*lay1.energy_bins_limits),
        lay1.zenith_bins_limits,
        lay1.detection_rate.T,
        norm=LogNorm(), vmin=1e-7, vmax=1e2   
    )
    plt.xlabel('log E/eV')
    plt.ylabel('Zenith [deg]')
    plt.title('triggered event rate over array '+'$\\nu_{ev}$'+', trihex all')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('$ [day^{-1} km^{-2}]$')
    plt.savefig(os.path.join(plot_path, "trig_ev_rate_trihexall.png"))

    plt.figure(110)
    plt.clf()
    plt.pcolor(
        np.log10(1e18*lay1.energy_bins_limits),
        lay1.zenith_bins_limits,
        lay1.differential_rate.T * 1e15,
        norm=LogNorm()  , vmin=1e-7, vmax=1e2
    )
    plt.xlabel('log E/eV')
    plt.ylabel('Zenith [deg]')
    plt.title('differential event rate over array '+'$\\nu_{ev}$'+', trihex all')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('$ [day^{-1} km^{-2} PeV^{-1} sr^{-1}]$')
    plt.savefig(os.path.join(plot_path, "diff_ev_rate_trihexall.png"))



    lay_tri1000.plot_2D_detection_rate()
    lay_tri1000.plot_2D_differential_rate()
