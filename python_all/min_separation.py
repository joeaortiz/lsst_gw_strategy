import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import random
#import os


from python_utilities.des_io import parse_observations
from parse_simlib import parse_simlib

libid_details, pointings = parse_simlib('../LSST_WFD_COADD.SIMLIB')
#libid_details islist with each row giving ['LIBID:', 'RA:', 'DECL:', 'NOBS:', 'MWEBV:', 'PIXSIZE:', 'first observation index', 'last observation index']


def find_ang_sep(RA1, DECL1, RA2, DECL2): #all inputs and output are in degrees and are floats
    event1 = SkyCoord(ra=float(RA1) * u.degree, dec=float(DECL1) * u.degree)
    event2 = SkyCoord(ra=float(RA2) * u.degree, dec=float(DECL2) * u.degree)
    angsep = event1.separation(event2)
    return angsep.deg


def get_timelocation(datfile):
    obs, headerdict = parse_observations(datfile)
    SIM_PEAKMJD = float(headerdict['SIM_PEAKMJD'])
    RA = float(headerdict['RA'])
    DECL = float(headerdict['DECL'])
    return SIM_PEAKMJD, RA, DECL


def day_min_separation(SIM_PEAKMJD, RA, DECL, band): #RA and DECL both in degrees
#returns the minimum angular separation and the details of the closest observation in the given band in the next 24 hours after the peakMJD
    relevant_pointings = [x for x in pointings if SIM_PEAKMJD < float(x[3]) and float(x[3]) < (SIM_PEAKMJD + 1) and x[5] == band] #pointings in the following 24hrs in the given band
    if relevant_pointings == []: #if no pointings in that band in the next 24hrs
        closest_pointing = []
        min_ang = np.nan
    else:
        #find the closest of observation
        ang_sep = np.zeros(len(relevant_pointings))
        for i, pointing in enumerate(relevant_pointings):
            ang_sep[i] = find_ang_sep(RA, DECL, pointing[1], pointing[2])
        argmin = np.argmin(ang_sep)
        min_ang = np.min(ang_sep)
        closest_pointing = relevant_pointings[argmin]
    return min_ang, closest_pointing



def make_sky_locations(num_locations):
    sky_locations = np.zeros([num_locations,3])
    for i in range(len(sky_locations)):
        sky_locations[i,:] = [random.uniform(59804,63180.9), random.uniform(0,360), random.uniform(-62,2)]
        #this range of values for DECL and PEAKMJD is justified in the notes in the googledoc
    return sky_locations




#DIR = '../GW170817_AT2017gfo_LSST_WFD/'
#numfiles = len([file for file in os.listdir(DIR) if file.endswith(".DAT")])

num_locations = 100
sky_locations = make_sky_locations(num_locations)
min_ang = np.zeros(num_locations)
bands = ['g', 'r', 'i', 'z', 'Y']
for band in bands:
    for i,location in enumerate(sky_locations):
        min_ang[i], closest_pointing = day_min_separation(location[0], location[1], location[2], band)
    print band, '\n', min_ang
    np.save('24hr_min_separation/min_ang_' + band, min_ang)