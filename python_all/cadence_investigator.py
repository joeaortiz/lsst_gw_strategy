import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

from python_utilities.des_io import parse_observations
from parse_simlib import parse_simlib


obs, headerdict = parse_observations('../GW170817_AT2017gfo_LSST_WFD/GW170817_AT2017gfo_LSST_WFD_SN000060.DAT')

libid_details, pointings = parse_simlib('../LSST_WFD_COADD.SIMLIB')
#libid_details islist with each row giving ['LIBID:', 'RA:', 'DECL:', 'NOBS:', 'MWEBV:', 'PIXSIZE:', 'first observation index', 'last observation index']
peak_MJD = float(headerdict['SIM_PEAKMJD'])

def find_ang_sep(RA1, DECL1, RA2, DECL2): #all inputs and output are in degrees and are floats
    event1 = SkyCoord(ra=float(RA1) * u.degree, dec=float(DECL1) * u.degree)
    event2 = SkyCoord(ra=float(RA2) * u.degree, dec=float(DECL2) * u.degree)
    angsep = event1.separation(event2)
    return angsep.deg


def find_closest(bin_pointings):
    dt_angsep = np.zeros([len(bin_pointings), 2])
    for i,pointing in enumerate(bin_pointings):
        dt_angsep[i,0] = float(pointing[1]) - peak_MJD
        dt_angsep[i,1] = find_ang_sep(float(headerdict['RA']), float(headerdict['DECL']), float(pointing[1]), float(pointing[2]))

    argmin = np.argmin(dt_angsep[:,1])
    DT_NEAR = dt_angsep[argmin, 0]
    ANGSEP_NEAR = dt_angsep[argmin, 1]
    if ANGSEP_NEAR != np.min(dt_angsep[:,1]):
        print 'error finding minimum angular separation'
    return DT_NEAR, ANGSEP_NEAR



binsize = 1/12. #1 day bins for now
nbins = 14*12
bands = ['g', 'r', 'i', 'z', 'Y']
plt.figure()
for band in bands:
    table = np.zeros([nbins, 3])
    for bin in range(nbins):
        MJD_start = peak_MJD + binsize*bin
        MJD_end = peak_MJD + binsize*(bin+1)
        bin_pointings = [x for x in pointings if MJD_start < float(x[3]) and float(x[3]) < MJD_end and x[5] == band]
        if len(bin_pointings) == 0:
            DT_NEAR, ANGSEP_NEAR = bin*binsize, 180
        else:
            DT_NEAR, ANGSEP_NEAR = find_closest(bin_pointings)
        table[bin] = [bin, DT_NEAR, ANGSEP_NEAR]
    print band
    print table
    table = table[table[:, 2] != 0., :]
    plt.plot(table[:,1], (180 - table[:,2]), marker='x', label=band)

plt.legend()
plt.xlabel('days after trigger')
plt.ylabel('180 - angular separation of nearest pointing')
plt.title('Angular separation of nearest pointing in 2 weeks after trigger, binsize=2hours')
plt.ylim(ymin=0)
plt.show()

