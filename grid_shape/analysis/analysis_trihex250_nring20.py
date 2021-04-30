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

path = "/Users/benoitl/Documents/GRAND/Data_grids/250trihex_nring20/"
plot_path = '/Users/benoitl/Documents/GRAND/Data_grids/250trihex_nring20/plots'
os.makedirs(plot_path, exist_ok=True)
primary = "Proton"
input_n_ring = 20
threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 10 # number of triggered antennas required to trigger an event
 
 
### "Creating" the layout
pos, offset = grids.create_grid_univ("trihex", 250, do_prune=False, input_n_ring=20)

## Creating the pruned layout
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 250, do_prune=True, input_n_ring=20)

# creating the all mask (all the antennas)
mask  = masks.make_all_mask(input_n_ring=20)

#print('toto')
if False:
    lay1 = layout.Layout(path, pos, mask, "all", threshold, n_trig_thres, input_n_ring, primary)
    lay1.make_all_plots()
    lay1.plot_2d_efficiency(6, 6, color='C0')

if False:
    lay2 = layout.Layout(path, pos2, mask2, "simple", threshold, n_trig_thres, input_n_ring, primary)
    lay2.make_all_plots()    
    lay2.plot_2d_efficiency(6, 6, color='C1')


if False:
    # creating the trihex n*250 grid out of the trixhex 250
    
    mask_tri500 = masks.make_coarse_out_of_fine_trihex(500, 'trihex', 250, 20) 
    lay_tri500 = layout.Layout(path, pos, mask_tri500, "tri500", threshold, n_trig_thres, input_n_ring, primary)
    lay_tri500.make_all_plots()
    lay_tri500.plot_2d_efficiency(6, 6, color='C2')

    mask_hex500 = masks.make_coarse_out_of_fine_trihex(500, 'hexhex', 250, 20) 
    lay_hex500 = layout.Layout(path, pos, mask_hex500, "hex500", threshold, n_trig_thres, input_n_ring, primary)
    lay_hex500.make_all_plots()
    lay_hex500.plot_2d_efficiency(6, 6, color='C3')

    mask_tri750 = masks.make_coarse_out_of_fine_trihex(750, 'trihex', 250, 20) 
    lay_tri750 = layout.Layout(path, pos, mask_tri750, "tri750", threshold, n_trig_thres, input_n_ring, primary)
    lay_tri750.make_all_plots()
    lay_tri750.plot_2d_efficiency(6, 6, color='C4')

    mask_hex750 = masks.make_coarse_out_of_fine_trihex(750, 'hexhex', 250, 20) 
    lay_hex750 = layout.Layout(path, pos, mask_hex750, "hex750", threshold, n_trig_thres, input_n_ring, primary)
    lay_hex750.make_all_plots()
    lay_hex750.plot_2d_efficiency(6, 6, color='C5')

    mask_tri1000 = masks.make_coarse_out_of_fine_trihex(1000, 'trihex', 250, 20) 
    lay_tri1000 = layout.Layout(path, pos, mask_tri1000, "tri1000", threshold, n_trig_thres, input_n_ring, primary)
    lay_tri1000.make_all_plots()
    lay_tri1000.plot_2d_efficiency(6, 6, color='C6')
if True:
    mask_hex1000 = masks.make_coarse_out_of_fine_trihex(1000, 'hexhex', 250, 20) 
    lay_hex1000 = layout.Layout(path, pos, mask_hex1000, "hex1000", threshold, n_trig_thres, input_n_ring, primary)
    lay_hex1000.make_all_plots()
    lay_hex1000.plot_2d_efficiency(6, 6, color='C7')

if False:
    mask_tri1250 = masks.make_coarse_out_of_fine_trihex(1250, 'trihex', 250, 20) 
    lay_tri1250 = layout.Layout(path, pos, mask_tri1250, "tri1250", threshold, n_trig_thres, input_n_ring, primary)
    lay_tri1250.make_all_plots()
    lay_tri1250.plot_2d_efficiency(6, 6, color='C6')

    mask_hex1250 = masks.make_coarse_out_of_fine_trihex(1250, 'hexhex', 250, 20) 
    lay_hex1250 = layout.Layout(path, pos, mask_hex1250, "hex1250", threshold, n_trig_thres, input_n_ring, primary)
    lay_hex1250.make_all_plots()
    lay_hex1250.plot_2d_efficiency(6, 6, color='C7')



if False:
    mask_1ring_250 = masks.make_coarse_out_of_fine_trihex_2(250, 'hexhex', 0, 125, 10)
    print('doing: lay_1ring_250')
    lay_1ring_250 = layout.Layout(
        path,
        pos,
        mask_1ring_250,
        "lay_1ring_250",
        threshold,
        n_trig_thres,
        input_n_ring=10,
        primary="Proton"
    )
    lay_1ring_250.make_all_plots()
    lay_1ring_250.plot_2d_efficiency(6, 6)



mask_trihex125_nring0 = masks.make_coarse_out_of_fine_trihex_2(125, 'trihex', 0, 125, 10)
mask_trihex125_nring1 = masks.make_coarse_out_of_fine_trihex_2(125, 'trihex', 1, 125, 10)
mask_trihex125_nring2 = masks.make_coarse_out_of_fine_trihex_2(125, 'trihex', 2, 125, 10)
mask_trihex125_nring3 = masks.make_coarse_out_of_fine_trihex_2(125, 'trihex', 3, 125, 10)
if False:
    print('doing: lay_trihex125_nring0')
    lay_trihex125_nring0 = layout.Layout(
            path,
            pos,
            mask_trihex125_nring0,
            "lay_trihex125_nring0",
            threshold,
            n_trig_thres,
            input_n_ring=10,
            primary="Proton"
        )
    lay_trihex125_nring0.make_all_plots()
    lay_trihex125_nring0.plot_2d_efficiency(6, 6, 'C3')

    print('doing: lay_trihex125_nring1')
    lay_trihex125_nring1 = layout.Layout(
            path,
            pos,
            mask_trihex125_nring1,
            "lay_trihex125_nring1",
            threshold,
            n_trig_thres,
            input_n_ring=10,
            primary="Proton"
        )
    lay_trihex125_nring1.make_all_plots()
    lay_trihex125_nring1.plot_2d_efficiency(6, 6, 'C4')
if False:
    print('doing: lay_trihex125_nring2')
    lay_trihex125_nring2 = layout.Layout(
            path,
            pos,
            mask_trihex125_nring2,
            "lay_trihex125_nring2",
            threshold,
            n_trig_thres,
            input_n_ring=10,
            primary="Proton"
        )
    lay_trihex125_nring2.make_all_plots()
    lay_trihex125_nring2.plot_2d_efficiency(6, 6, 'C5')

if False:
    print('doing: lay_trihex125_nring3')
    lay_trihex125_nring3 = layout.Layout(
            path,
            pos,
            mask_trihex125_nring3,
            "lay_trihex125_nring3",
            threshold,
            n_trig_thres,
            input_n_ring=10,
            primary="Proton"
        )
    lay_trihex125_nring3.make_all_plots()
    lay_trihex125_nring3.plot_2d_efficiency(6, 6, 'C6')

# lay_tri250.plot_layout()

# # creating the trihex 500 grid out of the trixhex 125
# mask_tri500 = masks.make_trihex_new_out_of_125(pos, 500, 10) 
# lay_tri500 = layout.Layout(path, pos, mask_tri500, "tri500", threshold, n_trig_thres, input_n_ring, primary)

# lay_tri500.plot_layout()



# # creating the trihex 500 grid out of the trixhex 125
# mask_tri1000 = masks.make_trihex_new_out_of_125(pos, 1000, 10) 
# lay_tri1000 = layout.Layout(path, pos, mask_tri1000, "tri1000", threshold, n_trig_thres, input_n_ring, primary)

# lay_tri1000.plot_layout()


# mask_island1 = masks.make_centralisland_out_of_125(pos, 10) 
# lay_island1 = layout.Layout(path, pos, mask_island1, "island1", threshold, n_trig_thres, input_n_ring, primary)

# lay_island1.plot_layout()



# mask_island2 = masks.make_centralisland_out_of_125_v2(pos, 10) 
# lay_island2 = layout.Layout(path, pos, mask_island2, "island2", threshold, n_trig_thres, input_n_ring, primary)

# lay_island2.plot_layout()



# mini_island_mask = masks.make_mini_island()
# lay_mini_island = layout.Layout(
#     path,
#     pos,
#     mini_island_mask,
#     "mini_island",
#     threshold, n_trig_thres, input_n_ring, primary
# )


# mini_island2_mask = masks.make_mini_island2()
# lay_mini_island2 = layout.Layout(
#     path,
#     pos,
#     mini_island2_mask,
#     "mini_island2",
#     threshold, n_trig_thres, input_n_ring, primary
# )






# plt.figure(3, figsize=(8,6))
# plt.clf()
# n_zenith_bins = len(lay1.zenith_bins_centers)
# n_bins_to_plot = 6
# i_bins = [(i) * n_zenith_bins // (n_bins_to_plot+1) for i in range(n_bins_to_plot)]
# #i_bins = [0, 1, 2, 3, 4]
# for k, i_bin in enumerate(i_bins):
#     plt.errorbar(
#         np.log10(lay1.energy_bins_centers*1e18),
#         lay1.detection_rate[:,i_bin],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=8,
#         color = "C1",
#         ls = '-',
#         label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[i_bin],lay1.zenith_bins_limits[i_bin+1])
#     )
#     plt.errorbar(
#         np.log10(lay1.energy_bins_centers*1e18),
#         lay_island1.detection_rate[:,i_bin],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=8,
#         color="C2",
#         ls = '-',
#         # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
#     )
#     plt.errorbar(
#         np.log10(lay1.energy_bins_centers*1e18),
#         lay_mini_island2.detection_rate[:,i_bin],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=8,
#         color="C3",
#         ls = '-',
#         # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
#     )
#     plt.errorbar(
#         np.log10(lay1.energy_bins_centers*1e18),
#         lay_tri500.detection_rate[:,i_bin],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=8,
#         color="C4",
#         ls = '-',
#         # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
#     )
#     plt.errorbar(
#         np.log10(lay1.energy_bins_centers*1e18),
#         lay_tri1000.detection_rate[:,i_bin],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=8,
#         color="C5",
#         ls = '-',
#         # label='%4.0f > zen >%4.0f deg'%(lay1.zenith_bins_limits[izen],lay1.zenith_bins_limits[izen+1])
#     )

#     plt.yscale('log')
#     plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} km^{-2}]$')
#     plt.xlabel('log E/eV')
#     legend1 = plt.legend(loc=1)
#     #splt.xscale('log')

# custom_lines = [
#     Line2D([0], [0], color="C1", lw=4),
#     Line2D([0], [0], color="C2", lw=4),
#     Line2D([0], [0], color="C3", lw=4),
#     Line2D([0], [0], color="C4", lw=4),
#     Line2D([0], [0], color="C5", lw=4)
# ]

# legend2 = plt.legend(custom_lines, ['trihex all', 'island 1',"mini island2", "trihex500", "trihex1000"], loc=3)

# plt.gca().add_artist(legend1)
# plt.gca().add_artist(legend2)
# plt.savefig(os.path.join(lay1.plot_path, "detection_rate_vs_energy.png"))




# plt.figure(4, figsize=(8,6))
# plt.clf()
# n_energy_bins = len(lay1.energy_bins_centers)
# n_bins_to_plot = 4
# i_bins = [(i+1) * n_energy_bins // (n_bins_to_plot+1) for i in range(n_bins_to_plot)]
# #i_bins = [0, 1, 2, 3, 4]
# for k, i_bin in enumerate(i_bins):
# #for iener in range(0, len(lay1.energy_bins_limits)-1):
#     plt.errorbar(
#         lay1.zenith_bins_centers,
#         lay1.detection_rate[i_bin,:],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=10,
#         ls = '-',
#         color="C1",
#         label='%4.2f > log E/eV >%4.2f'%(
#             np.log10(lay1.energy_bins_limits[i_bin]*1e18),
#             np.log10(lay1.energy_bins_limits[i_bin+1]*1e18)
#         )
#     )
#     plt.errorbar(
#         lay1.zenith_bins_centers,
#         lay_island1.detection_rate[i_bin,:],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=10,
#         ls = '-',
#         color="C2"
#         )
#     plt.errorbar(
#         lay_mini_island2.zenith_bins_centers,
#         lay_mini_island2.detection_rate[i_bin,:],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=10,
#         ls = '-',
#         color="C3"
#     )
#     plt.errorbar(
#         lay_tri500.zenith_bins_centers,
#         lay_tri500.detection_rate[i_bin,:],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=10,
#         ls = '-',
#         color="C4"
#     )
#     plt.errorbar(
#         lay_tri1000.zenith_bins_centers,
#         lay_tri1000.detection_rate[i_bin,:],
#         fmt=SYM_LIST[i_bin],
#         capsize=2,
#         alpha=0.7,
#         ms=10,
#         ls = '-',
#         color="C5"
#     )


#     plt.yscale('log')
#     plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} km^{-2}]$')
#     plt.xlabel('Zenith [deg]')
#     #plt.title('%s, %4.2f > zenith > %4.2f deg'%(lay1.mask_name, lay1.[izen], thetal[izen+1]))
#     legend1 = plt.legend(loc=1)
    
 
# custom_lines = [
#     Line2D([0], [0], color="C1", lw=4),
#     Line2D([0], [0], color="C2", lw=4),
#     Line2D([0], [0], color="C3", lw=4),
#     Line2D([0], [0], color="C4", lw=4),
#     Line2D([0], [0], color="C5", lw=4)
# ]

# legend2 = plt.legend(custom_lines, ['trihex all', 'island1',"mini island2", "trihex500", "trihex1000"], loc=3)

# plt.gca().add_artist(legend1)
# plt.gca().add_artist(legend2)
# plt.savefig(os.path.join(lay1.plot_path, "detection_rate_vs_zenith.png"))


# plt.figure(45)
# lay1.plot_layout()
# plt.savefig(os.path.join(lay1.plot_path, "lay1.png"))



# plt.figure(46)
# lay_rand_5.plot_layout(fig=46)
# plt.savefig(os.path.join(lay1.plot_path, "lay_rand_5.png"))


# plt.figure(47)
# lay_tri250.plot_layout(fig=47)
# plt.savefig(os.path.join(lay1.plot_path, "lay_tri250.png"))


# plt.figure(48)
# lay_tri500.plot_layout(fig=48)
# plt.savefig(os.path.join(lay1.plot_path, "lay_tri500.png"))



# plt.figure(49)
# lay_tri1000.plot_layout(fig=49)
# plt.savefig(os.path.join(lay1.plot_path, "lay_tri1000.png"))


# plt.figure(50)
# lay_island1.plot_layout(fig=50)
# plt.savefig(os.path.join(lay1.plot_path, "lay_island1.png"))


# plt.figure(51)
# lay_island2.plot_layout(fig=51)
# plt.savefig(os.path.join(lay1.plot_path, "lay_island2.png"))


# plt.figure(52)
# lay_mini_island.plot_layout(fig=51)
# plt.savefig(os.path.join(lay1.plot_path, "lay_mini_island.png"))


# plt.figure(53)
# lay_mini_island2.plot_layout(fig=53)
# plt.savefig(os.path.join(lay1.plot_path, "lay_mini_island2.png"))




# lay1.make_diff_event_rate(200)
# lay_tri1000.make_diff_event_rate(200)
# lay_tri500.make_diff_event_rate(200)
# lay_tri250.make_diff_event_rate(200)
# lay_island1.make_diff_event_rate(200)

# lay_mini_island.make_diff_event_rate(200)
# lay_mini_island2.make_diff_event_rate(200)


# plt.figure(345)
# plt.clf()
# plt.plot(lay1.energy_bins_centers, lay1.integrated_ev_rate_no_trig, label="trig efficiency = 1")
# plt.plot(lay1.energy_bins_centers, lay1.integrated_ev_rate, label="trihex all")
# plt.plot(lay_tri250.energy_bins_centers, lay_tri250.integrated_ev_rate, label="trihex250")
# plt.plot(lay_tri500.energy_bins_centers, lay_tri500.integrated_ev_rate, label="trihex500")
# plt.plot(lay_tri1000.energy_bins_centers, lay_tri1000.integrated_ev_rate, label="trihex1000")
# plt.plot(lay_island1.energy_bins_centers, lay_island1.integrated_ev_rate, label="island1")
# plt.plot(lay_mini_island.energy_bins_centers, lay_mini_island.integrated_ev_rate, label="mini island")
# plt.plot(lay_mini_island2.energy_bins_centers, lay_mini_island2.integrated_ev_rate, label="mini island2")

# plt.ylabel('Differential event rate [day'+"$^{-1}$"+"PeV"+"$^{-1}$"+']')
# plt.title('Detector area 200 km'+"$^2$")
# plt.xlabel('Energy [EeV]')
# plt.xscale('log')
# plt.yscale('log')
# plt.legend(loc=0)
# plt.savefig('diff_rate_200km2.png')




# ####   hist 2d




# lay1.plot_2D_detection_rate()
# lay1.plot_2D_differential_rate()
# lay1.plot_mean_n_trig()


# lay_tri250.plot_2D_detection_rate()
# lay_tri250.plot_2D_differential_rate()
# lay_tri250.plot_mean_n_trig()


# lay_tri500.plot_2D_detection_rate()
# lay_tri500.plot_2D_differential_rate()
# lay_tri500.plot_mean_n_trig()


# lay_tri1000.plot_2D_detection_rate()
# lay_tri1000.plot_2D_differential_rate()
# lay_tri1000.plot_mean_n_trig()


# lay_island1.plot_2D_detection_rate()
# lay_island1.plot_2D_differential_rate()
# lay_island1.plot_mean_n_trig()


# lay_island2.plot_2D_detection_rate()
# lay_island2.plot_2D_differential_rate()
# lay_island2.plot_mean_n_trig()


# lay_mini_island2.plot_2D_detection_rate()
# lay_mini_island2.plot_2D_differential_rate()
# lay_mini_island2.plot_mean_n_trig()

