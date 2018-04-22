import subprocess
import os
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import scipy
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/data/des41.a/data/jaortiz/snana-test/')
from python_all.parse_simlib import parse_simlib
os.chdir('/data/des41.a/data/jaortiz/snana-test/do_swaps/')
from parse_dump import findall_KN 



def amake_swaps(t1, t2, max_airmass, max_extra_slew, pointings_arr):
    os.chdir('/data/des41.a/data/jaortiz/snana-test/')
    print t1, t2, max_airmass, max_extra_slew
    
    success = []
    maybe_obs = []
    
    plt.figure()

    
    allKN = findall_KN()
    KNlist = allKN
    
    os.chdir('/data/des41.a/data/jaortiz/snana-test/')

    #Creates a list of all pointings, with [LIBID, RA, DECL, MJD IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG]

    #Results analysis 
    #[swap made, observation already scheduled, no observations in next 24hrs, slew angle >70deg, no observations in next week close to KN]

    bands = ['u', 'g', 'r', 'i', 'z', 'Y']
    new_pointings_arr = np.copy(pointings_arr)

    slews = np.zeros([len(KNlist),6])
#     no_obs, obs_made = np.zeros(6, dtype=float), np.zeros(6, dtype=float)
    costs = [ [], [], [], [], [], [] ]

    folder = str(t1) + '_' + str(t2) + '_' + str(max_airmass) + '_' + str(max_extra_slew)

    subprocess.call(["mkdir", "do_swaps/new_pointings_arrays/" + folder])

    
    #logging
    os.chdir('/data/des41.a/data/jaortiz/snana-test/do_swaps/new_pointings_arrays/' + folder + '/')
    log = open('swaps.log', 'w')
    log.write('t1 = %i, t2 = %i, max airmass = %i, max slew angle = %i \n\n' %(t1, t2, max_airmass, max_extra_slew))
    log.write('Number of KN processed: %i \n\n' %len(KNlist))
    log.write('Swap made, KN already scheduled to be observed in next t1 days, KN is not visible in next t1 days, No observations when KN is visible in next t1 days, Closest observation has max slew angle too large, No observation in next t2 days close to KN \n')

    results = np.zeros([6,8]) #must reset variable for each band

    for j, KN in enumerate(KNlist):
        #print '*****************New KN****************', KN
        if j%100 ==0:
            print j
        
        success_count = 0
        
        pointings_later = pointings_arr[KN[0]<pointings_arr[:,3], :]
        pointings_nextday = pointings_later[pointings_later[:,3]<(KN[0]+t1)]
        #print 'Number of pointings in next t1 days: ', len(pointings_nextday)
        y1 = np.full(len(pointings_nextday), 7)
        plt.scatter(pointings_nextday[:,3], y1, marker=1)
        
        obs_times = when_observable(KN, max_airmass, t1)
        #print 'Number of observable times: ', len(obs_times)
        y = np.zeros(len(obs_times))
        plt.scatter(obs_times, y, marker='x')

        for band in range(6):
            #print '%s Filter' %bands[band]
            band_pointings_arr = pointings_arr[pointings_arr[:,5]==band, :]

            #when is the KN going to be visible? Check the next 24 hours in 0.5 hour intervals


            #Are we scheduled to observe this event anyway? If not find the closest observation in the next 24 hours
            closest_obs, slew_to_KN, original_overhead_time = get_closest_obs(pointings_arr, band_pointings_arr, KN, obs_times, t1, max_extra_slew, band)
            slews[j, band] = slew_to_KN

            if type(original_overhead_time) == int:
                if original_overhead_time == 0:
                    results[band,3] += 1 #no observations in next 24hrs
                    #print 'Swap failed as no observations when KN is visible in next %i day(s)' %t1
                    continue
                if original_overhead_time == 1:
                    results[band,1] += 1 #observation already exists
                    #print 'No swap required as observation already exists with slew to KN: %f' %slew_to_KN
                    success_count += 1
                    continue
                if original_overhead_time == 2:
                    results[band,4] += 1 #slew angle > max_extra_slew deg
                    #print 'Swap failed as closest observation has angsep: %f > 70deg' %slew_to_KN
                    continue
                if original_overhead_time == 3:
                    results[band,2] += 1 #slew angle > max_extra_slew deg
                    #print 'Swap failed as KN not visible in next %i day(s)' %t1
                    continue

            #Search for future pointing at KN location (within 0.5deg) in next 7 days.
            #if we have made it this far then we want to do a swap unless there is no observation in the future we can swap with
            pointing_to_swap, min_angsep_lt = find_later_pointing(band_pointings_arr, KN, t2)
            if type(pointing_to_swap) == int:
                results[band,5] += 1 #note that no swap because there is no pointing close to KN in next week
                #print 'Swap failed as no observation in next %i days close to KN. Closest obs at %f' %(t2, min_angsep_lt)
                continue

            #What is the cost in slew time of swapping these pointings?
            extra_overhead_time, swap_ang = cost(pointing_to_swap, closest_obs, original_overhead_time, pointings_arr, j, KN)
            costs[band].append(extra_overhead_time)

            #swap libid, ra and decl 
            new_pointings_arr = swap_obs(new_pointings_arr, closest_obs, pointing_to_swap)
            results[band,0] += 1
            #print j, 'Swap done. t1swap= ', closest_obs[3]-KN[0], 'swap_angle= ', swap_ang
            success_count += 1

        plt.yticks(np.arange(8),['KN is visible', 'u', 'g', 'r', 'i', 'z', 'Y', 'observation in any band'])
        #plt.show()
        
        if success_count ==1:
            maybe_obs.append(KN[0])
        if success_count > 1:
            success.append(KN[0])

    results_frac = results/float(len(KNlist))
#     no_obs[band] = results_frac[2] + results_frac[3] #to display this figure on graph too
#     obs_made[band] = results_frac[1] 
    for band in range(len(bands)):
        log.write('%s: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f \n' %(bands[band],results_frac[band,0],results_frac[band,1],results_frac[band,2],\
                                                                results_frac[band,3],results_frac[band,4],results_frac[band,5]) )

    #        print '\n', '****************** %s Filter *********************' %bands[band]
    #        print 'Number of KN processed: %i \n' %len(KNlist)
    #        print 'Swap made: %.1f  \nNo swap made because... \n    Observation already scheduled: %.1f ' %(results_frac[0], results_frac[1])
    #        print '    No observations in next %i day(s): %.1f ' %(t1,results_frac[3])
    #        print '    Swap failed as extra slew required > %i deg: %.1f ' %(max_extra_slew,results_frac[4])
    #        print '    No observation in next %i days close to KN: %.1f ' %(t2, results_frac[5])
    #        print  '***************************************************'

    #     log.close()

    print 'Num observed once: ', len(maybe_obs)
    print 'Num KN observed twice or more: ', len(success)
        
    data = [slews[:,0][~np.isnan(slews[:,0])], slews[:,1][~np.isnan(slews[:,1])], slews[:,2][~np.isnan(slews[:,2])], \
                 slews[:,3][~np.isnan(slews[:,3])], slews[:,4][~np.isnan(slews[:,4])], slews[:,5][~np.isnan(slews[:,5])] ]
    

    
    
    tot = 0
    for c in costs:
        c = np.array(c)
        #print len(cost), sum(np.array(cost)<0), sum(cost)
        pos = c[c>0]*24*3600
        neg = c[c<0]*24*3600
        tot += sum(pos)
    
    summary = np.array([t1, t2, max_airmass, max_extra_slew, 0,0,0,0,tot, len(success)])

    os.chdir('/data/des41.a/data/jaortiz/snana-test/')
    np.save("do_swaps/new_pointings_arrays/" + folder + '/new_pointings_arr', new_pointings_arr)
    np.save("do_swaps/new_pointings_arrays/" + folder + '/costs', costs)
    np.save("do_swaps/new_pointings_arrays/" + folder + '/success', success)
    np.save("do_swaps/new_pointings_arrays/" + folder + '/data', data)
    np.save("do_swaps/new_pointings_arrays/" + folder + '/summary', summary)

    return new_pointings_arr, data, costs, success


def when_observable(KN, max_airmass, t1):#when is the KN going to be visible? Check the next 24 hours in 0.5 hour intervals

    lsst_location = EarthLocation(lat=-30.2446*u.deg, lon=-70.7494*u.deg, height=2663*u.m)
    times = np.arange(KN[0], KN[0]+t1, 1/48.) #30 minute intervals in next 24 hours
    time = Time(times, format='mjd')
    KN_loc = SkyCoord(ra=KN[1] * u.degree, dec=KN[2] * u.degree)
    KNaltaz = KN_loc.transform_to(AltAz(obstime=time,location=lsst_location))  
    altitude = KNaltaz.alt.deg
    
    #put limit on airmass, Could be below ~2
    zenith_ang = 90 - altitude #in deg
    zenith_ang_rad = zenith_ang * np.pi / 180
    airmass = 1 / np.cos(zenith_ang_rad)
    
    observable = (airmass<max_airmass).astype(int) #observable if altitude in sky is >5deg
    for i,t in enumerate(times):
        if airmass[i] <0:
            observable[i]=0
    
    #ind = np.nonzero(times)
    obs_times = observable*times
    ind = np.array(observable.astype(bool))
    obs_times = obs_times[ind]
    
    return obs_times 


def get_closest_obs(pointings_arr, band_pointings_arr, KN, obs_times, t1, max_extra_slew, band):
    #Are we scheduled to observe this event anyway? If not find the closest observation in the next 24 hours

    #get pointings in the next week
    pointings_later = band_pointings_arr[KN[0]<band_pointings_arr[:,3], :]
    pointings_nextday_possible = pointings_later[pointings_later[:,3]<(KN[0]+t1)]
    
    #print '%i possible pointings in next t1 days' %len(pointings_nextday_possible)#, pointings_nextday_possible[:,3]
    if len(pointings_nextday_possible) != 0:
        y = np.full(len(pointings_nextday_possible[:,3]), band+1)
        plt.scatter(pointings_nextday_possible[:,3], y, marker='x')
    
    #only consider pointings in next day when KN is observable
    if len(obs_times) == 0: #if there are no observable times then we still need to create variable pointings_nextday
        pointings_nextday = []
        closest_obs, original_overhead_time, slew_to_KN = 0,3, np.nan
    else:
        for i,time in enumerate(obs_times):
            add_maybe = pointings_nextday_possible[pointings_nextday_possible[:,3]<(time+ 1/48.)]
            add = add_maybe[add_maybe[:,3]>time]
            if i==0:
                pointings_nextday = np.copy(add)
            else:
                pointings_nextday = np.vstack((pointings_nextday, add))
    
        #pointings_not_observable = [pointing for pointing in pointings_nextday_possible if pointing not in pointings_nextday]
        #print len(pointings_nextday_possible), len(pointings_not_observable), len(pointings_nextday)
    if len(pointings_nextday_possible) != 0: #for checking    
        pointings_lcs = SkyCoord(ra=pointings_nextday_possible[:,1] * u.degree, dec=pointings_nextday_possible[:,2] * u.degree)
        KN_loc = SkyCoord(ra=KN[1] * u.degree, dec=KN[2] * u.degree)
        ags = KN_loc.separation(pointings_lcs).deg
                
    
    if len(pointings_nextday) == 0:
        if len(obs_times) != 0 : #if no pointings in the next 24hrs, due to no pointings
            closest_obs, original_overhead_time, slew_to_KN = 0,0, np.nan
    else:
        
        #print 'Number of pointings in next day when KN is observable', len(pointings_nextday[:,3])

        #get angles between all of pointings and KN
        pointings_locs = SkyCoord(ra=pointings_nextday[:,1] * u.degree, dec=pointings_nextday[:,2] * u.degree)
        KN_loc = SkyCoord(ra=KN[1] * u.degree, dec=KN[2] * u.degree)
        angsep = KN_loc.separation(pointings_locs).deg
        
        #find closest pointing within 24hours for swap
        closest_obs = pointings_nextday[np.argmin(angsep),:]
        slew_to_KN = np.min(angsep)
        
        row_of_pointingsarr = np.argwhere(pointings_arr[:,4]==closest_obs[4])[0][0]
        original_overhead_time = pointings_arr[row_of_pointingsarr,3] - pointings_arr[row_of_pointingsarr-1,3]
        
        #check if we are scheduled to observe this event anyway
        #observations_arg = np.argwhere(angsep < np.sqrt(9.6/np.pi))
        #KN_observations = pointings_nextday[observations_arg,:] #observations made of KN, if any
        if slew_to_KN < np.sqrt(9.6/np.pi):
            original_overhead_time = 1 #no swap because observation already exists
        if slew_to_KN > max_extra_slew:
            original_overhead_time = 2 #no swap because smallest slew angle was too large

    if len(pointings_nextday_possible) != 0:
        if np.min(ags) < np.sqrt(9.6/np.pi) and np.min(ags) != slew_to_KN:
            print 'Problem', np.min(ags), slew_to_KN
            mjd = pointings_nextday_possible[np.argmin(ags), 3]
            print mjd, obs_times
            lsst_location = EarthLocation(lat=-32.344633333333334*u.deg, lon=-77.34941666666666*u.deg, height=2652*u.m)
            times = np.arange(KN[0], KN[0]+t1, 1/48.) #30 minute intervals in next 24 hours
            time = Time(mjd, format='mjd')
            KN_loc = SkyCoord(ra=KN[1] * u.degree, dec=KN[2] * u.degree)
            KNaltaz = KN_loc.transform_to(AltAz(obstime=time,location=lsst_location))  
            altitude = KNaltaz.alt.deg

            #put limit on airmass, Could be below ~2
            zenith_ang = 90 - altitude #in deg
            zenith_ang_rad = zenith_ang * np.pi / 180
            airmass = 1 / np.cos(zenith_ang_rad)
            print 'airmass is ', airmass
    
    return closest_obs, slew_to_KN, original_overhead_time



def find_later_pointing(pointings_arr, KN, t2):
    #Search for future pointing at KN location (within 0.5deg) in next 7 days.
    #if we have made it this far then we want to do a swap unless there is no observation in the future we can swap with

    #get pointings in the next week
    pointings_later = pointings_arr[KN[0]<pointings_arr[:,3], :]
    pointings_nextweek = pointings_later[pointings_later[:,3]<(KN[0]+t2)]
    
    nextweek_locs = SkyCoord(ra=pointings_nextweek[:,1] * u.degree, dec=pointings_nextweek[:,2] * u.degree)
    KN_loc = SkyCoord(ra=KN[1] * u.degree, dec=KN[2] * u.degree)
    angsep_longterm = KN_loc.separation(nextweek_locs).deg
    min_angsep_lt = np.min(angsep_longterm)
    if np.min(angsep_longterm) > np.sqrt(9.6/np.pi): #Require KN to be inside observation as we will swap ra,decl
        pointing_to_swap = 0 #no observations in next week to swap with
    else:
        possible_swap_arg = np.argwhere(angsep_longterm < np.sqrt(9.6/np.pi)) 
        pointing_to_swap = pointings_nextweek[possible_swap_arg,:][0][0] #get pointing to swap
        
    return pointing_to_swap, min_angsep_lt




def cost(pointing_to_swap, closest_obs, original_overhead_time, pointings_arr, j, KN):
    #What is the cost in slew time of swapping these pointings?

    swap_to_now = SkyCoord(ra=closest_obs[1] * u.degree, dec=closest_obs[2] * u.degree)
    swap_to_later = SkyCoord(ra=pointing_to_swap[1] * u.degree, dec=pointing_to_swap[2] * u.degree)
    swap_angle = swap_to_now.separation(swap_to_later).deg 
    const1, const2 = 33.5250434/(24*3600.), 1.67705303/(24*3600.) #these numbers must be changed
    new_overhead_time = const1 + const2*swap_angle
    extra_overhead_time = new_overhead_time - original_overhead_time #units of days
    if original_overhead_time <0:
        print 'original overhead time is less than zero!'
    if extra_overhead_time < -0.00025:
        print 'extra overhead time is negative and large:', extra_overhead_time, 'New: ', new_overhead_time*3600*24, 'Old: ', original_overhead_time*3600*24
        row_of_pointingsarr = np.argwhere(pointings_arr[:,4]==closest_obs[4])[0][0]
        print j, closest_obs[3], pointings_arr[row_of_pointingsarr-2:row_of_pointingsarr+3,3]
    #if extra_overhead_time > 0.0005:
    #    print 'extra overhead time is large. New: ', new_overhead_time*3600*24, 'Old: ', original_overhead_time*3600*24

        
    return extra_overhead_time, swap_angle


def swap_obs(new_pointings_arr, closest_obs, pointing_to_swap): #swap libid, ra and decl 
    
    for i in range(len(new_pointings_arr)):
        if new_pointings_arr[i,4]==closest_obs[4]:
            new_pointings_arr[i,0:3] = pointing_to_swap[0:3]
        if new_pointings_arr[i,4]==pointing_to_swap[4]:
            new_pointings_arr[i,0:3] = closest_obs[0:3]

    return new_pointings_arr
