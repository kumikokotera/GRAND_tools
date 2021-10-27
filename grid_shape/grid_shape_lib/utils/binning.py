import nunmpy as np


def compute_energy_bins(
    low_energy=16.3,
    high_energy=18.7,
    delta_log10=0.1
):

    dum_e = np.arange(low_energy, high_energy+delta_log10/10, delta_log10) 
    energy_bins_limits = 10**(
        dum_e - delta_log10/2
    ) / 1e18

    energy_bins_centers = 10**(dum_e[0:-1])/1e18
    delta_energy = energy_bins_limits[1:] - energy_bins_limits[:-1]
    return energy_bins_limits, energy_bins_centers, delta_energy


def compute_zenith_bins(
    sec_theta_min=1.162,
    sec_theta_max=21.54,
    n_zenith_bins=16
):
    log_sec_theta_min = np.log10(sec_theta_min)
    log_sec_theta_max = np.log10(sec_theta_max)
    log_increment = (log_sec_theta_max-log_sec_theta_min)/n_zenith_bins
    # print("On log spacing")
    # print("max log_sec_theta:"+str(max_theta)+" -> " + str(log_sec_theta_max))
    # print("min log_sec_theta:"+str(min_theta)+" -> " + str(log_sec_theta_min))
    # print("log increment is sec theta:"+str(log_increment)+" a factor " + str(np.power(10,log_increment)))
    limits_log_secant_bins = np.logspace(
        log_sec_theta_min,
        log_sec_theta_max,
        n_zenith_bins+1,
        base=10
    )
    # print("log_limits:" + str(limits_log_secant_bins))
    center_log_secant_bins = limits_log_secant_bins[0:n_zenith_bins]*np.power(10,log_increment/2)
    # print("log centers:"+ str(center_log_secant_bins))
    zenith_bins_limits = np.rad2deg(np.arccos(1.0/limits_log_secant_bins))
    # print("limits:" + str(self.zenith_bins_limits))
    zenith_bins_centers = np.rad2deg(np.arccos(1.0/center_log_secant_bins))
    # print("centers:" + str(self.zenith_bins_centers))
    return zenith_bins_limits, zenith_bins_centers
