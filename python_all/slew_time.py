import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import scipy

#Creates a list of all pointings, with [LIBID, RA, DECL, MJD IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG]
from parse_simlib import parse_simlib

libid_details, pointings_arr = parse_simlib('../minion_1016_WFD.simlib')


#slew_times = [overhead time, angular separation, change band?] changeband = 0 for no change and 1 for change of band.
#length of overhead_times is len(mjd_band) - 1
slew_times = np.zeros([len(pointings_arr)-1, 3])
for i in range(len(slew_times)):
    slew_times[i,0] = np.round(pointings_arr[i+1,3] - pointings_arr[i,3], 4)
    if pointings_arr[i+1,5] != pointings_arr[i,5]:
        slew_times[i,2] = 1
    
    #get angular separation of the pointings to the KN
    loc1 = SkyCoord(ra=pointings_arr[i,1] * u.degree, dec=pointings_arr[i,2] * u.degree)
    loc2 = SkyCoord(ra=pointings_arr[i+1,1] * u.degree, dec=pointings_arr[i+1,2] * u.degree)
    slew_times[i,1] = loc1.separation(loc2).deg
    
    if i%10000 == 0:
        print i
    
np.save('slew_data', slew_times)