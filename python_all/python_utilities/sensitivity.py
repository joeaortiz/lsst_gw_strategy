from collections import OrderedDict as odict

import cosmology
import cosmology_constants
import matplotlib
import numpy as np
import pylab
import scipy.interpolate
import star_formation_rate

import des_io
import snana_fitsio

pylab.ion()

############################################################

params = {
    #'backend': 'eps',
    'axes.labelsize': 16,
    #'text.fontsize': 12,           
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    'text.usetex': True,
    #'figure.figsize': fig_size,
    'font.family':'serif',
    'font.serif':'Computer Modern Roman',
    'font.size': 12
    }
matplotlib.rcParams.update(params)

############################################################

def detectEvents(listfile, obs):
    cid_array, detect_array = snana_fitsio.detectAll(listfile)

    master_detect_array = np.tile(False, len(obs))
    select = np.in1d(obs['CID'], cid_array[detect_array])
    master_detect_array[select] =  True
    detect_array = master_detect_array
    
    return detect_array

############################################################

def detectionEfficiency(obs, detect_array, plot=False):
    zbins = np.arange(0., 1. + 1.e-10, 0.05)
    efficiency_array = np.zeros(len(zbins) - 1)
    efficiency_weight_array = np.zeros(len(zbins) - 1)
    zcenter_array = 0.5 * (zbins[0:-1] + zbins[1:])

    cut_base = (obs['PEAKMJD'] >= 56550) & (obs['PEAKMJD'] <= 56579) & (obs['GENTYPE'] != 1)
    for ii in range(0, len(zbins) - 1):
        zmin = zbins[ii]
        zmax = zbins[ii + 1]
    
        cut_z = (obs['GENZ'] > zmin) & (obs['GENZ'] < zmax)

        n_total = np.sum(cut_base & cut_z)
        n_detect = np.sum(cut_base & cut_z & detect_array)
    

        weight_total = np.sum((cut_base & cut_z) * 10**(-0.4 * obs['MAGT0_r']))
        weight_detect = np.sum((cut_base & cut_z & detect_array) * 10**(-0.4 * obs['MAGT0_r']))
    
        print zmin, zmax, n_total, n_detect
    
        if n_total > 0:
            efficiency_array[ii] = float(n_detect) / n_total
            efficiency_weight_array[ii] = weight_detect / weight_total
        else:
            efficiency_array[ii] = -1.
            efficiency_weight_array[ii] = -1.

    cut = (efficiency_array >= 0.)
    zcenter_array_clean = np.concatenate([[0.], zcenter_array[cut], [2.]])
    efficiency_array_clean = np.concatenate([[1.], efficiency_array[cut], [0.]])
    f_efficiency = scipy.interpolate.interp1d(zcenter_array_clean, efficiency_array_clean, bounds_error=False, fill_value=0.)

    cut = (efficiency_weight_array >= 0.)
    zcenter_weight_array_clean = np.concatenate([[0.], zcenter_array[cut], [2.]])
    efficiency_weight_array_clean = np.concatenate([[1.], efficiency_weight_array[cut], [0.]])
    f_efficiency_weight = scipy.interpolate.interp1d(zcenter_weight_array_clean, efficiency_weight_array_clean, bounds_error=False, fill_value=0.)
    
    if plot:
        z = np.linspace(0., 1., 1000)
        
        pylab.figure()
        pylab.plot(z, f_efficiency(z), c='red', lw=2, label='Uniform Weight')
        pylab.plot(z, f_efficiency_weight(z), c='green', lw=2, label='Weighted by Peak Optical Luminosity')
        pylab.scatter(zcenter_array, efficiency_array, c='red', edgecolor='none')
        pylab.scatter(zcenter_array, efficiency_weight_array, c='green', edgecolor='none')
        pylab.ylim(0., 1.)
        pylab.xlim(0., 1.)
        pylab.xlabel('Redshift')
        pylab.ylabel('Detection Efficiency')
        pylab.legend(loc='upper right')
        if save:
            pylab.savefig('detection_efficiency.pdf')

    return f_efficiency, f_efficiency_weight

############################################################

def cumulativeIntensity(f_efficiency, f_efficiency_weight, plot=False):
    z_array = np.linspace(1.e-4, 10., 10000)
    dz = z_array[1] - z_array[0]

    # SN rate (s^-1 sr^-1)
    rate_array = cosmology.D_H * cosmology.D_C(z_array)**2 \
                 * (1. + z_array)**(-1) * cosmology.E(z_array)**(-1) \
                 * star_formation_rate.snrMadau(z_array) \
                 * dz

    # With time dilation
    intensity_array = (cosmology_constants.C_LIGHT / (cosmology_constants.HUBBLE * 4 * np.pi)) \
                      * (1. + z_array)**(-3) * cosmology.E(z_array)**(-1) \
                      * star_formation_rate.snrMadau(z_array) \
                      * dz

    detected_intensity_array = (cosmology_constants.C_LIGHT / (cosmology_constants.HUBBLE * 4 * np.pi)) \
                               * (1. + z_array)**(-3) * cosmology.E(z_array)**(-1) \
                               * star_formation_rate.snrMadau(z_array) \
                               * f_efficiency(z_array) \
                               * dz

    detected_weight_intensity_array = (cosmology_constants.C_LIGHT / (cosmology_constants.HUBBLE * 4 * np.pi)) \
                                      * (1. + z_array)**(-3) * cosmology.E(z_array)**(-1) \
                                      * star_formation_rate.snrMadau(z_array) \
                                      * f_efficiency_weight(z_array) \
                                      * dz

    if plot:
        pylab.figure()
        #pylab.plot(z_array, intensity_array)
        pylab.plot(z_array, np.cumsum(intensity_array) / np.sum(intensity_array), lw=2, c='blue', label='Total Emission')
        pylab.plot(z_array, np.cumsum(detected_intensity_array) / np.sum(intensity_array), lw=2, c='red', label='Detected Emission (uniform)')
        pylab.plot(z_array, np.cumsum(detected_weight_intensity_array) / np.sum(intensity_array), lw=2, c='green', label='Detected Emission (weighted)')
        pylab.xlabel('z')
        pylab.ylabel('Neutrino Intensity CDF')
        pylab.xlim(0., 4.)
        pylab.ylim(0., 1.)
        pylab.legend(loc='lower right')
        if save:
            pylab.savefig('neutrino_intensity_cdf.pdf', dpi=200)

        pylab.xlim(0., 1.)
        pylab.ylim(0., 0.5)
        if save:
            pylab.savefig('neutrino_intensity_cdf_zoom.pdf', dpi=200)

    return z_array, \
        np.cumsum(intensity_array) / np.sum(intensity_array), \
        np.cumsum(detected_intensity_array) / np.sum(intensity_array), \
        np.cumsum(detected_weight_intensity_array) / np.sum(intensity_array)

############################################################

def backgroundRate(obs, detect_array):
    # Estimate background rate

    AREA_SIM = 300.
    #area_icecube_conservative = np.pi * 1.0**2
    area_icecube_conservative = 3.
    area_icecube_optimistic = np.pi * 0.5**2

    cut_dict = odict()
    cut_dict['30_all'] = (obs['PEAKMJD'] >= 56550) & (obs['PEAKMJD'] <= 56579)
    cut_dict['10_all'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569)
    cut_dict['10_cc'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENTYPE'] != 1)
    #cut_dict['10_photoz'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENZ'] < 0.5)
    #cut_dict['10_photoz_cc'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENZ'] < 0.5) & (obs['GENTYPE'] != 1)
    cut_dict['10_photoz'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENZ'] < 0.4)
    cut_dict['10_photoz_cc'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENZ'] < 0.4) & (obs['GENTYPE'] != 1)
    cut_dict['10_specz_cc'] = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENZ'] < 0.2) & (obs['GENTYPE'] != 1)

    for key in cut_dict.keys():
        print 'Cut = %s'%(key)
        print '  Rate = %.2f'%(np.sum(detect_array[cut_dict[key]]) / AREA_SIM)
        print '  Background = %.2f'%((area_icecube_conservative / AREA_SIM) * np.sum(detect_array[cut_dict[key]]))
        print '  Background = %.2f'%((area_icecube_optimistic / AREA_SIM) * np.sum(detect_array[cut_dict[key]]))

    z = np.sort(obs['GENZ'][detect_array & cut_dict['10_photoz_cc']])
    rate = np.linspace(0.,
                       (area_icecube_conservative / AREA_SIM) * np.sum(detect_array[cut_dict['10_photoz_cc']]),
                       np.sum(detect_array[cut_dict['10_photoz_cc']]))
    #pylab.figure()
    #pylab.plot(z, r)
    return z, rate

############################################################

def detectedEvents(obs, detect_array):
    cut_ia_detected = (obs['GENTYPE'] == 1) & detect_array
    cut_cc_detected = (obs['GENTYPE'] != 1) & detect_array
    pylab.figure()
    pylab.scatter(obs['GENZ'][cut_ia_detected], obs['MAGT0_r'][cut_ia_detected], c='black', marker='o', s=5, edgecolor='none')
    pylab.scatter(obs['GENZ'][cut_cc_detected], obs['MAGT0_r'][cut_cc_detected], c='red', marker='o', s=10, edgecolor='none')
    pylab.xlim(0., 1.)
    pylab.ylim(16., 25.)
    pylab.xlabel('Redshift')
    pylab.ylabel('Peak r-band Magnitude')
    
    print 'MEDIAN =', np.median(obs['MAGT0_r'][cut_cc_detected & (obs['GENZ'] < 0.4)])

############################################################

def eventRateSummary(obs, detect_array, f_efficiency, z_array, cum_intensity):
    
    AREA_SIM = 300.
    AREA_DECAM_POINTING = 3.
    N_FOLLOWUP = 4. # 4 per semester
    F_PURITY =  0.5

    cut_signal = (obs['PEAKMJD'] >= 56550) & (obs['PEAKMJD'] <= 56579) & (obs['GENTYPE'] != 1)
    print 'CUT_SIGNAL'
    print np.sum(cut_signal)
    
    # TRYING TO COME UP WITH A MAGNITUDE DISTRIBUTION OF SIGNAL EVENTS
    #n_signal = 1000
    #f_redshift_inverse = scipy.interpolate.interp1d(cum_intensity, z_array, bounds_error=False, fill_value=0.)
    #z_randoms = f_redshift_inverse(np.random.random(n_signal))
    #obs['GENZ']

    cut_10_all = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & detect_array
    cut_10_photoz_cc = (obs['PEAKMJD'] >= 56560) & (obs['PEAKMJD'] <= 56569) & (obs['GENZ'] < 0.5) & (obs['GENTYPE'] != 1) & detect_array

    print 'CUT_10_ALL'
    print np.sum(cut_10_all)

    weights_10_all = np.tile(N_FOLLOWUP * AREA_DECAM_POINTING / AREA_SIM, 
                             np.sum(cut_10_all))
    weights_10_photoz_cc = np.tile(N_FOLLOWUP * AREA_DECAM_POINTING / AREA_SIM, 
                             np.sum(cut_10_photoz_cc))

    bins = np.arange(16., 25.1, 0.25)
    pylab.figure()
    pylab.hist(obs['MAGT0_r'][cut_10_all], bins=bins, color='0.5', histtype='step', lw=2, cumulative=True, weights=weights_10_all)
    pylab.hist(obs['MAGT0_r'][cut_10_photoz_cc], bins=bins, color='black', histtype='step', lw=2, cumulative=True, weights=weights_10_photoz_cc)
    pylab.xlim(16., 25.1)
    pylab.xlabel('Peak $r$-band Magnitude')
    pylab.ylabel('Cumulative Number of Events Per DES Season')

############################################################


save = False

results = {}
for maglim in [21, 24]:
    datadir = '/Users/keithbechtol/Documents/DES/projects/icecube/snana/v2_fits/KB_DES+ICECUBE_SN_m%2i_FITS'%(maglim)
    dumpfile = '%s/KB_DES+ICECUBE_SN_m%2i_FITS.DUMP'%(datadir, maglim)
    listfile = '%s/KB_DES+ICECUBE_SN_m%2i_FITS.LIST'%(datadir, maglim)
    print datadir
    
    params, obs = des_io.parse_dump(dumpfile)

    detect_array = detectEvents(listfile, obs)
    f_efficiency, f_efficiency_weight = detectionEfficiency(obs, detect_array)
    z_array, cum_intensity, cum_detected_intensity, cum_detected_intensity_weight = cumulativeIntensity(f_efficiency, f_efficiency_weight)

    results['%2i'%(maglim)] = {'z_array': z_array,
                               'cum_intensity': cum_intensity, 
                               'cum_detected_intensity': cum_detected_intensity,
                               'cum_detected_intensity_weight': cum_detected_intensity_weight}

    z_background, rate_background = backgroundRate(obs, detect_array)
    results['%2i'%(maglim)]['z_background'] = z_background
    results['%2i'%(maglim)]['rate_background'] = rate_background

    detectedEvents(obs, detect_array)
    #eventRateSummary(obs, detect_array)

pylab.figure()
total, = pylab.plot(results['21']['z_array'], results['21']['cum_intensity'], lw=2, c='black', label='From All SN')
#pylab.plot(np.NaN, np.NaN, '-', color='none', label='w/ Optically Detected Counterparts')
first_legend = pylab.legend(handles=[total], loc='upper left')
pylab.gca().add_artist(first_legend)
unc_21, = pylab.plot(results['21']['z_array'], results['21']['cum_detected_intensity'], lw=2, c='0.5', label=r'$m_{lim}$ = 21, uncorrelated')
wht_21, = pylab.plot(results['21']['z_array'], results['21']['cum_detected_intensity_weight'], lw=2, c='0.5', ls='--', label=r'$m_{lim}$ = 21, weighted')
unc_24, = pylab.plot(results['24']['z_array'], results['24']['cum_detected_intensity'], lw=2, c='red', label=r'$m_{lim}$ = 24, uncorrelated')
wht_24, = pylab.plot(results['24']['z_array'], results['24']['cum_detected_intensity_weight'], lw=2, c='red', ls='--', label=r'$m_{lim}$ = 24, weighted')
pylab.xlabel('Redshift')
#pylab.ylabel('Cumulative Neutrino Intensity')
pylab.ylabel('Cumulative Fraction of Total Neutrino Intensity')
#pylab.xlim(0., 4.)
#pylab.ylim(0., 1.)
pylab.xlim(0., 1.)
pylab.ylim(0., 0.5)
pylab.legend(handles=[unc_21, wht_21, unc_24, wht_24], loc='lower right', title='From Optically Detectable SN')
if True:
    pylab.savefig('neutrino_intensity_cdf_21_24_compare_analysis.pdf')
    pylab.savefig('neutrino_intensity_cdf_21_24_compare_analysis.eps')

#pylab.xlim(0., 1.)
#pylab.ylim(0., 0.5)
#if save:
#    pylab.savefig('neutrino_intensity_cdf_zoom.pdf', dpi=200)

#####

N_EVENTS = 4.
PURITY = 0.5
pylab.figure()
total, = pylab.plot(results['21']['z_array'], N_EVENTS * PURITY * results['21']['cum_intensity'], lw=2, c='black', label='Total Associated')
#first_legend = pylab.legend(handles=[total], loc='upper left')
#pylab.gca().add_artist(first_legend)
unc_21, = pylab.plot(results['21']['z_array'], N_EVENTS * PURITY * results['21']['cum_detected_intensity'], lw=2, c='0.5', label=r'Associated, $m_{lim}$ = 21')
unc_24, = pylab.plot(results['24']['z_array'], N_EVENTS * PURITY * results['24']['cum_detected_intensity'], lw=2, c='red', label=r'Associated, $m_{lim}$ = 24')
br_24, = pylab.plot(results['24']['z_background'], N_EVENTS * results['24']['rate_background'], lw=2, c='blue', ls='--', label=r'Unassociated, $m_{lim}$ = 24')
pylab.xlabel('Redshift')
pylab.ylabel('Cumulative Number of Optically Detected CC SN')
pylab.title('Follow-up for 4 IceCube Alerts in 2017B\n50\% Purity of Astrophysical Neutrino Sample')
#pylab.xlim(0., 4.)
#pylab.ylim(0., 1.)
pylab.xlim(0., 0.4)
#pylab.ylim(0., 1.5)
pylab.ylim(0., 1.0)
#pylab.legend(loc='upper left', frameon=False)
pylab.legend(loc='upper right', frameon=False)
if save:
    pylab.savefig('detection_cdf_21_24_compare_analysis.pdf')
    pylab.savefig('detection_cdf_21_24_compare_analysis.eps')

