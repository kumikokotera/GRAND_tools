import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import numpy as np
import os
import json
from grid_shape_lib.utils import diff_spec as diff_spec
from grid_shape_lib.modules import grids as grids
from grid_shape_lib.utils import utils_analysis as ua
from grid_shape_lib.modules import layout as layout
from grid_shape_lib.modules import masks as masks


import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)
matplotlib.rc('legend', fontsize=18)

matplotlib.rcParams['font.family'] = "sans-serif"

SYM_LIST = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
MYC = ['0','0.20','0.4','0.6','0.8']


### Commented analysis script


path = "/Users/ab212678/Documents/GRAND/Data_grids/3_dec_20"
plot_path = '/Users/ab212678/Documents/GRAND/Data_grids/3_dec_20/plots_GP100_jan23_plot_for_paper'


os.makedirs(plot_path, exist_ok=True)
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/Trihex"

## use 75 !!!!
threshold = 75  # trigger threshold for individual antennas in muV
n_trig_thres = 5  # number of triggered antennas required to trigger an event
 
 
### "Creating" the layout
pos, offset = grids.create_grid_univ("trihex", 250, do_prune=False, input_n_ring=20)
## Creating the pruned layout
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 250, do_prune=True, input_n_ring=20)


# creating the all mask (all the antennas)
mask = masks.make_all_mask(input_n_ring=20)



mask_tri1250 = masks.make_coarse_out_of_fine_trihex(1250, "trihex", 250, 20)
mask_tri1000 = masks.make_coarse_out_of_fine_trihex(1000, "trihex", 250, 20)
mask_tri500 = masks.make_coarse_out_of_fine_trihex(500, "trihex", 250, 20)
mask_tri250 = masks.make_coarse_out_of_fine_trihex(250, "trihex", 250, 20)
mask_tri750 = masks.make_coarse_out_of_fine_trihex(750, "trihex", 250, 20)

mask_hex1250 = masks.make_coarse_out_of_fine_trihex(1250, "hexhex", 250, 20)
mask_hex1000 = masks.make_coarse_out_of_fine_trihex(1000, "hexhex", 250, 20)
mask_hex750 = masks.make_coarse_out_of_fine_trihex(750, "hexhex", 250, 20)
mask_hex500 = masks.make_coarse_out_of_fine_trihex(500, "hexhex", 250, 20)
mask_hex250 = masks.make_coarse_out_of_fine_trihex(250, "hexhex", 250, 20)


lay_all = layout.Layout(path, pos, mask, "all", threshold, n_trig_thres, input_n_ring=20, primary="Proton", plot_path=plot_path)
#lay_all.make_all_plots(do_title=False)



lay_tri1250 = layout.Layout(
    path,
    pos,
    mask_tri1250,
    "lay_tri1250",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_tri1250.make_all_plots(do_title=False)



lay_tri1000 = layout.Layout(
    path,
    pos,
    mask_tri1000,
    "lay_tri1000",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_tri1000.make_all_plots(do_title=False)


lay_tri500 = layout.Layout(
    path,
    pos,
    mask_tri500,
    "lay_tri500",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_tri500.make_all_plots(do_title=False)

lay_tri750 = layout.Layout(
    path,
    pos,
    mask_tri750,
    "lay_tri750",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_tri750.make_all_plots(do_title=False)

lay_tri250 = layout.Layout(
    path,
    pos,
    mask_tri250,
    "lay_tri250",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_tri250.make_all_plots(do_title=False)




lay_hex1250 = layout.Layout(
    path,
    pos,
    mask_hex1250,
    "lay_hex1250",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_hex1250.make_all_plots(do_title=False)



lay_hex1000 = layout.Layout(
    path,
    pos,
    mask_hex1000,
    "lay_hex1000",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_hex1000.make_all_plots(do_title=False)

lay_hex500 = layout.Layout(
    path,
    pos,
    mask_hex500,
    "lay_hex500",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_hex500.make_all_plots(do_title=False)


lay_hex750 = layout.Layout(
    path,
    pos,
    mask_hex750,
    "lay_hex750",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_hex750.make_all_plots(do_title=False)


lay_hex250 = layout.Layout(
    path,
    pos,
    mask_hex250,
    "lay_hex250",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_hex250.make_all_plots(do_title=False)







mask100_step1250 = masks.make_coarse_out_of_fine_trihex_2(1250, 'hexhex', 3, 250, 20)
masks.plot_mask(mask100_step1250, 20, "mask100_step1500")
lay_mask100_step1250 = layout.Layout(
    path,
    pos,
    mask100_step1250,
    "mask100_step1250",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_mask100_step1250.make_all_plots(do_title=False)


mask100_step1250 = masks.make_coarse_out_of_fine_trihex_2(1250, 'hexhex', 3, 250, 20)
R1250_3= (1 + 3 * 1.5) * 2 / np.sqrt(3) * 1250
lay_mask100_step1250_contained = layout.Layout(
    path,
    pos,
    mask100_step1250,
    "mask100_step1250_contained",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path,
    radius_contained=R1250_3
)
#lay_mask100_step1250_contained.make_all_plots(do_title=False)


mask100_step1000 = masks.make_coarse_out_of_fine_trihex_2(1000, 'hexhex', 3, 250, 20)
masks.plot_mask(mask100_step1000, 20, "mask100_step1000")
lay_mask100_step1000 = layout.Layout(
    path,
    pos,
    mask100_step1000,
    "mask100_step1000",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_mask100_step1000.make_all_plots(do_title=False)

R1000_3 = (1 + 3 * 1.5) * 2 / np.sqrt(3) * 1000
lay_mask100_step1000_contained = layout.Layout(
    path,
    pos,
    mask100_step1000,
    "mask100_step1000_contained",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path,
    radius_contained=R1000_3
)
#lay_mask100_step1000_contained.make_all_plots(do_title=False)

mask100_step750 = masks.make_coarse_out_of_fine_trihex_2(750, 'hexhex', 3, 250, 20)
#masks.plot_mask(mask100_step750, 20, "mask100_step750")
lay_mask100_step750 = layout.Layout(
    path,
    pos,
    mask100_step750,
    "mask100_step750",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_mask100_step750.make_all_plots(do_title=False)

R750_3 = (1 + 3 * 1.5) * 2 / np.sqrt(3) * 750
lay_mask100_step750_contained = layout.Layout(
    path,
    pos,
    mask100_step750,
    "mask100_step750_contained",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path,
    radius_contained=R750_3
)
#lay_mask100_step750_contained.make_all_plots(do_title=False)


mask100_step500 = masks.make_coarse_out_of_fine_trihex_2(500, 'hexhex', 3, 250, 20)
#masks.plot_mask(mask100_step500, 20, "mask100_step500")
lay_mask100_step500 = layout.Layout(
    path,
    pos,
    mask100_step500,
    "mask100_step500",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_mask100_step500.make_all_plots(do_title=False)

R500_3 = (1 + 3 * 1.5) * 2 / np.sqrt(3) * 500
lay_mask100_step500_contained = layout.Layout(
    path,
    pos,
    mask100_step500,
    "mask100_step500_contained",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path,
    radius_contained=R500_3
)
#lay_mask100_step500_contained.make_all_plots(do_title=False)

mask100_step250 = masks.make_coarse_out_of_fine_trihex_2(250, 'hexhex', 3, 250, 20)
#masks.plot_mask(mask100_step250, 20, "mask100_step250")
lay_mask100_step250 = layout.Layout(
    path,
    pos,
    mask100_step250,
    "mask100_step250",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
#lay_mask100_step250.make_all_plots(do_title=False)

R250_3 = (1 + 3 * 1.5) * 2 / np.sqrt(3) * 250
lay_mask100_step250_contained = layout.Layout(
    path,
    pos,
    mask100_step250,
    "mask100_step250_contained",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path,
    radius_contained=R250_3
)
#lay_mask100_step250_contained.make_all_plots(do_title=False)








list_to_plot = [
    lay_hex1000,
    lay_tri1000,
    lay_hex500,
    lay_tri500
]
labels = [
    "$\lambda$ = 1000 m, hex",
    "$\lambda$ = 1000 m, tri",
    "$\lambda$ = 500 m, hex",
    "$\lambda$ = 500 m, tri",
]

fmts = [
    'k-',
    'k:',
    'r-', 
    'r:'
]

# list_to_plot = [
#     # lay_mask100_step1250,
#     # lay_mask100_step1250_contained,
#     # lay_mask100_step1000,
#     # lay_mask100_step1000_contained,
#     lay_mask100_step750,
#     lay_mask100_step750_contained,
#     lay_mask100_step500,
#     lay_mask100_step500_contained,
#     # lay_mask100_step250,
#     # lay_mask100_step250_contained,

# ]






# plt.figure(figsize=(13, 10))
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step250_contained.detection_rate_summed_vs_zenith/lay_mask100_step250.detection_rate_summed_vs_zenith, 'C0', label='250m'
# )
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step250_contained.area/lay_mask100_step250.area*lay_all.zenith_bins_centers/lay_all.zenith_bins_centers, 'C0--', label='250m'
# )

# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step500_contained.detection_rate_summed_vs_zenith/lay_mask100_step500.detection_rate_summed_vs_zenith, 'C1', label='500m'
# )
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step500_contained.area/lay_mask100_step500.area*lay_all.zenith_bins_centers/lay_all.zenith_bins_centers, 'C1--', label='500m'
# )

# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step750_contained.detection_rate_summed_vs_zenith/lay_mask100_step750.detection_rate_summed_vs_zenith, 'C2', label='750m'
# )
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step750_contained.area/lay_mask100_step750.area*lay_all.zenith_bins_centers/lay_all.zenith_bins_centers, 'C2--', label='500m'
# )

# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step1000_contained.detection_rate_summed_vs_zenith/lay_mask100_step1000.detection_rate_summed_vs_zenith, 'C3', label='1000m'
# )
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step1000_contained.area/lay_mask100_step1000.area*lay_all.zenith_bins_centers/lay_all.zenith_bins_centers, 'C3--', label='500m'
# )
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step1250_contained.detection_rate_summed_vs_zenith/lay_mask100_step1250.detection_rate_summed_vs_zenith, 'C4',  label='1250m'
# )
# plt.plot(
#     lay_all.zenith_bins_centers,
#     lay_mask100_step1250_contained.area/lay_mask100_step1250.area*lay_all.zenith_bins_centers/lay_all.zenith_bins_centers, 'C4--', label='500m'
# )
# plt.legend(loc=0)
# plt.ylabel('ratio of rates and areas')
# plt.xlabel('Zenith [deg]')
# plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
# plt.xscale('log')

# plt.legend(loc=0, ncol=2)
# plt.tight_layout()
# plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_ratio_containedGP100_jan23.png'))








plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.zenith_bins_centers, lay.detection_rate_summed_vs_zenith, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=3, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_1000_500_tri_hex.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.energy_bins_centers, lay.detection_rate_summed_vs_energy, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_1000_500_tri_hex.png'))






plt.figure(figsize=(13, 10))
for lay in list_to_plot:
    plt.plot(lay.energy_bins_centers, lay.detection_rate_summed_vs_energy, label=lay.mask_name)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]")
plt.xlabel('Energy [EeV]')
plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0, ncol=2)
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_GP100_jan23.png'))


plt.figure(34545, figsize=(8, 6))
plt.clf()
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_energy, lay.energy_bins_limits, label=lab, color=fmt[0], ls=fmt[1:], lw=3)
#     plt.bar(lay.energy_bins_centers,lay.detection_rate_summed_vs_energy, lay.delta_energy, fill=False, log=True, label=lay.mask_name)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_1000_500_tri_hex_step.png'))


plt.figure(3454445, figsize=(8, 6))
plt.clf()
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_zenith, lay.zenith_bins_limits, label=lab, color=fmt[0], ls=fmt[1:], lw=3)
#     plt.bar(lay.energy_bins_centers,lay.detection_rate_summed_vs_energy, lay.delta_energy, fill=False, log=True, label=lay.mask_name)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=3, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_1000_500_tri_hex_step.png'))



plt.figure(43567)
plt.clf()
#plt.stairs(lay_hex1000.detection_rate_summed_vs_energy/lay_hex500.detection_rate_summed_vs_energy, lay_hex1000.energy_bins_limits, label='hex1000/hex500', lw=2)
#plt.stairs(lay_tri1000.detection_rate_summed_vs_energy/lay_tri500.detection_rate_summed_vs_energy, lay_tri1000.energy_bins_limits, label='tri1000/tri500', lw=2)
plt.stairs(lay_tri1000.detection_rate_summed_vs_energy/lay_hex1000.detection_rate_summed_vs_energy, lay_tri1000.energy_bins_limits,  label='$\lambda=1000$ m', color='k',  lw=2)
plt.stairs(lay_tri500.detection_rate_summed_vs_energy/lay_hex500.detection_rate_summed_vs_energy, lay_tri1000.energy_bins_limits, label='$\lambda=500$ m', lw=2, color='r')

plt.xlabel('Energy [EeV]')
plt.ylabel('Ratio $N_{ev}^{tri}/N_{ev}^{hex}$')
plt.legend(loc=0, ncol=1)
plt.xscale('log')
plt.ylim(1, 5)
plt.xlim(5e-2, 5)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'ratio_detection_rate_vs_energy_v3.png'))


############# plots de l'impact du step size ###########






# if False:

#     mask_mini_islands = masks.make_mini_islands(250, 20)
#     mask_1ring_250 = masks.make_coarse_out_of_fine_trihex_2(250, 'hexhex', 0, 250, 20)
#     mask_1ring_500 = masks.make_coarse_out_of_fine_trihex_2(500, 'hexhex', 0, 250, 20)
#     mask_1ring_1000 = masks.make_coarse_out_of_fine_trihex_2(1000, 'hexhex', 0, 250, 20)



   

#     mask_hex1250 = masks.make_coarse_out_of_fine_trihex(1250, "hexhex", 250, 20)







list_to_plot = [
    lay_hex1250,
    #lay_tri1250,
    lay_hex1000,
    #lay_tri1000,
    lay_hex750,
    #lay_tri750,
    lay_hex500,
    #lay_tri500,
    lay_hex250,
    #lay_tri250
]

labels = [
    "$\lambda$ = 1250 m, hex",
    #"$\lambda$ = 1250 m, tri",
    "$\lambda$ = 1000 m, hex",
    #"$\lambda$ = 1000 m, tri",
    "$\lambda$ = 750 m, hex",
    #"$\lambda$ = 750 m, tri",
    "$\lambda$ = 500 m, hex",
    #"$\lambda$ = 500 m, tri",
    "$\lambda$ = 250 m, hex",
    #"$\lambda$ = 250 m, tri",
]

fmts = [
    'k-',
    #'k:',
    'r-', 
    #'r:',
    'g-',
    #'g:',
    'b-', 
    #'b:',
    'm-', 
    #'m:',
]




plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.zenith_bins_centers, lay.detection_rate_summed_vs_zenith, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=3, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_1250_1000_750_500_250_tri_hex.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.energy_bins_centers, lay.detection_rate_summed_vs_energy, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_1250_1000_750_500_250_tri_hex.png'))





plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_zenith, lay.zenith_bins_limits,  color=fmt[0], ls=fmt[1:], lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=3, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_1250_1000_750_500_250_tri_hex_stairs.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_energy, lay.energy_bins_limits, color=fmt[0], ls=fmt[1:],  lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_1250_1000_750_500_250_tri_hex_stairs.png'))


























lay_hex1000_5ant_75muv = layout.Layout(
    path,
    pos,
    mask_hex1000,
    "lay_hex1000",
    75,
    5,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_hex1000_5ant_75muv.make_all_plots(do_title=False)

lay_hex1000_5ant_150muv = layout.Layout(
    path,
    pos,
    mask_hex1000,
    "lay_hex1000",
    150,
    5,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_hex1000_5ant_150muv.make_all_plots(do_title=False)

lay_hex1000_10ant_75muv = layout.Layout(
    path,
    pos,
    mask_hex1000,
    "lay_hex1000",
    75,
    10,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_hex1000_10ant_75muv.make_all_plots(do_title=False)

lay_hex1000_10ant_150muv = layout.Layout(
    path,
    pos,
    mask_hex1000,
    "lay_hex1000",
    150,
    10,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_hex1000_10ant_150muv.make_all_plots(do_title=False)






list_to_plot = [
    lay_hex1000_5ant_75muv,
    lay_hex1000_5ant_150muv,
    lay_hex1000_10ant_75muv,
    lay_hex1000_10ant_150muv
]

labels = [
    "$N_{\\rm{trig}}$ = 5, $V_{\\rm{trig}}$ = 75 $\mu$V",
    "$N_{\\rm{trig}}$ = 5, $V_{\\rm{trig}}$ = 150 $\mu$V",
    "$N_{\\rm{trig}}$ = 10, $V_{\\rm{trig}}$ = 75 $\mu$V",
    "$N_{\\rm{trig}}$ = 10, $V_{\\rm{trig}}$ = 150 $\mu$V",
]
# labels = [
#     "$N_{ant}$ = 5, thr= 75$\mu$V, $\lambda$ = 1000 m, hex",
#     "$N_{ant}$ = 5, thr= 150$\mu$V, $\lambda$ = 1000 m, hex",
#     "$N_{ant}$ = 10, thr= 75$\mu$V, $\lambda$ = 1000 m, hex",
#     "$N_{ant}$ = 10, thr= 150$\mu$V, $\lambda$ = 1000 m, hex",
# ]

fmts = [
    'b-',
    #'k:',
    'b--', 
    #'r:',
    'm-',
    #'g:',
    'm--', 
]



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.zenith_bins_centers, lay.detection_rate_summed_vs_zenith, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$'+', number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_ant_thr_hex.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.energy_bins_centers, lay.detection_rate_summed_vs_energy, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_ant_thr_hex.png'))







plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_zenith, lay.zenith_bins_limits, color=fmt[0], ls=fmt[1:], lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$'+', number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=2, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_ant_thr_hex_stairs.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs( lay.detection_rate_summed_vs_energy, lay.energy_bins_limits, color=fmt[0], ls=fmt[1:], lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=3, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_ant_thr_hex_stairs.png'))
























mask_hex1000_island250 =  masks.make_central_island_out_of_fine(
    1000, "hexhex", 250, "hexhex", 3, 250, 20
)
lay_hex1000_island250 = layout.Layout(
        path,
        pos,
        mask_hex1000_island250,
        "lay_hex1000_island250",
        threshold,
        n_trig_thres,
        input_n_ring=20,
        primary="Proton",
        plot_path=plot_path
    )
lay_hex1000_island250.make_all_plots(do_title=False)


mask_limited_mini_island = masks.make_limited_mini_island(
    'hexhex',
    250,
    1000,'hexhex',"trihex", 6, 250, 250, 20)

lay_limited_mini_island = layout.Layout(
    path,
    pos,
    mask_limited_mini_island,
    "lay_limited_mini_island",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_limited_mini_island.make_all_plots(do_title=False)


mask_spiral7 = masks.make_spiral_mask_with_minipose(1000, 1.8, 65, "hexhex", minipos_nring=0, minipos_shape='trihex')
#masks.make_spiral_mask_with_minipose(1000, 2.3, 20, "hexhex", minipos_nring=1, minipos_shape='hexhex')
masks.plot_mask(mask_spiral7, 20, "spiral7")
lay_spiral7 = layout.Layout(
    path,
    pos,
    mask_spiral7,
    "spiral7",
    threshold,
    n_trig_thres,
    input_n_ring=20,
    primary="Proton",
    plot_path=plot_path
)
lay_spiral7.make_all_plots(do_title=False)






list_to_plot = [
    lay_tri1000,
    lay_hex1000,
    lay_hex1000_island250,
    lay_limited_mini_island,
    lay_spiral7,
   
]

labels = [
    "$\lambda$ = 1000 m, tri",
    "$\lambda$ = 1000 m, hex",
    "island",
    "flower",
    "spiral",
    
]

fmts = [
    'r-',
    'k-',
    'g-', 
    #'r:',
    'b-',
    #'g:',
    'm-', 
]

### line plots 

plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.zenith_bins_centers, lay.detection_rate_summed_vs_zenith, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$'+', number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=4, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_infill_withhex1000.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.plot(lay.energy_bins_centers, lay.detection_rate_summed_vs_energy, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=1)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_infill_withhex1000.png'))


### staire plots


plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_zenith, lay.zenith_bins_limits, label=lab, color=fmt[0], ls=fmt[1:], lw=3)
    #plt.plot(lay.zenith_bins_centers, lay.detection_rate_summed_vs_zenith, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$'+', number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Zenith [deg]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
#plt.xscale('log')
plt.yscale('log')
plt.xlim(40, 90)
plt.legend(loc=4, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_zenith_paper_infill_withhex1000_stairs.png'))



plt.figure(figsize=(8, 6))
for lay, lab, fmt in zip(list_to_plot, labels, fmts):
    plt.stairs(lay.detection_rate_summed_vs_energy, lay.energy_bins_limits, label=lab, color=fmt[0], ls=fmt[1:], lw=3)
    #plt.plot(lay.energy_bins_centers, lay.detection_rate_summed_vs_energy, fmt, lw=3, label=lab)
plt.ylabel('$N_{\\rm{ev}}$, '+'number of events [day'+"$^{-1}$]", fontsize=20)
plt.xlabel('Energy [EeV]', fontsize=20)
#plt.title('Detection rate, {} muV, {} antennas'.format(threshold, n_trig_thres))
plt.xscale('log')
plt.yscale('log')
plt.xlim(5e-2, 5)
plt.legend(loc=0, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(plot_path, 'day_detection_rate_vs_energy_paper_infill_withhex1000_stairs.png'))




### saving stuff for Xin Xu

np.savetxt('lay_hex1000_zeniths_bin_centers.txt', lay_hex1000.zenith_bins_centers)
np.savetxt('lay_hex1000_detection_rate_summed_vs_zenith.txt', lay_hex1000.detection_rate_summed_vs_zenith)