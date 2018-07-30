"""This investigates the time taken for a given KN event to be observed."""


import numpy as np
import random
from astropy import units as u
from astropy.coordinates import SkyCoord
import scipy
import os
import matplotlib.pylab as plt

from parse_dump import findall_KN

os.chdir('/data/des41.a/data/jaortiz/snana-test/python_all/')
from parse_simlib import parse_simlib
from python_utilities.des_io import parse_observations

def make_histogram_vars(simlib, all_KN): 

    #Creates list of all pointings with [LIBID, RA, DECL, MJD IDEXPT FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG]
    libid_details, pointings = parse_simlib(simlib)
    print 'Number of LSST pointings: ', len(pointings)
    
    if all_KN == 'no':
    
        DIR = '/data/des41.b/data/SNDATA_ROOT/SIM/GW170817_AT2017gfo_LSST_WFD/'
        #get KN locations of simulated events in snana
        file_count = len([f for f in os.walk(DIR).next()[2] if f[-4:] == ".DAT"])
        KN_locations = np.zeros([file_count, 3])

        # run through all .dat files in the specified directory and get time and location of KN
        file_num = 0
        for file in os.listdir(DIR):  
            if file.endswith(".DAT"):
                obs, headerdict = parse_observations(DIR + file)
                KN_locations[file_num,:] = [headerdict['PEAKMJD'], headerdict['RA'], headerdict['DECL']]
                file_num += 1    
        print 'Number of KN observed by snana: ', len(KN_locations)
    
    elif all_KN == 'yes':
        
        KN_locations = findall_KN()
     
    else:
        print 'Invalid choice entry for all_KN'
    
    
    #get time taken until observed for all KN
    time_taken = np.zeros(len(KN_locations))
    time_taken[:] = np.nan
    for i, KN in enumerate(KN_locations):
        if i%50 == 0:
            print i

        #get pointings in the next week
        pointings_later = pointings[KN[0]<pointings[:,3], :]
        pointings_nextweek = pointings_later[pointings_later[:,3]<(KN[0]+7)]
        if len(pointings_nextweek) == 0: #if no pointings in the next 24hrs
            continue

        #get angles between all of pointings and KN
        kilonova_loc = SkyCoord(ra=KN[1] * u.degree, dec=KN[2] * u.degree)
        pointings_locs = SkyCoord(ra=pointings_nextweek[:,1] * u.degree, dec=pointings_nextweek[:,2] * u.degree)
        angsep = kilonova_loc.separation(pointings_locs).deg
        observations_arg = np.argwhere(angsep<np.sqrt(9.6/np.pi))
        if len(observations_arg) == 0:
            continue
        mjd_observed = pointings_nextweek[observations_arg[0],:][0][3]

        time_taken[i] = mjd_observed - KN[0]
        
    return time_taken, KN_locations
    
    
def plot_histogram(time_taken, KN_locations, all_KN, x_max):
    
    efficiency = len(KN_locations)/11659.
    
    no_observation_nextweek = len(time_taken[np.isnan(time_taken)])
    fraction_no_obs_nextweek = float(no_observation_nextweek)/ len(KN_locations) 
    print 'fraction_no_obs_nextweek: ', fraction_no_obs_nextweek
       
    #plot histogram
    histdata, binedges = np.histogram(time_taken[~np.isnan(time_taken)], 25)
    binsize = binedges[1] - binedges[0]
    bincentres = binedges[:-1] + binsize/2
    histdata = histdata/float(len(time_taken[~np.isnan(time_taken)])) 

    fig1 = plt.figure(figsize=(12,10))
    plt.subplot(2,1,1)
    plt.plot(bincentres, np.cumsum(histdata)*efficiency*(1-fraction_no_obs_nextweek), label='cumulative histogram')
    plt.plot(bincentres, histdata*efficiency*(1-fraction_no_obs_nextweek), alpha=0.7, label='histogram')
    if all_KN == 'no':
        plt.title('Histogram of time taken until first observation of Kilonova. Efficiency = %f' % efficiency, size=16, fontweight='bold')
    elif all_KN == 'yes':
        plt.title('Histogram of time taken until first observation of Kilonova.', size=16, fontweight='bold')
    plt.ylabel('Fraction of KN', size=12)
    plt.legend(loc= 'center right')
    plt.xlim(xmax=x_max)


    plt.subplot(2,1,2)
    plt.plot(bincentres, histdata*efficiency*(1-fraction_no_obs_nextweek), alpha=0.7, label='histogram', color='g')
    plt.ylabel('Fraction of KN', size=12)
    plt.xlabel('Time after KN of first observation (days)', size=12)
    plt.xlim(xmax=x_max)


    #plt.hist(time_taken[~np.isnan(time_taken)], bins=40)
    plt.show()
    
    return 