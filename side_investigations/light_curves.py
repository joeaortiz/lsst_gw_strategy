"""Make light curve plot for a given .dat file"""



import numpy as np
import matplotlib.pyplot as plt
import os


from python_utilities.des_io import parse_observations


def lightcurves(datfile): #makes lightcurve plot for a given .dat file
    obs, headerdict = parse_observations(datfile)

    bands = [] #get list of bands
    for row in obs:
        if row[1] not in bands:
            bands.append(row[1])

    plt.figure(datfile)
    for band in bands:
        MJD, FLUXCAL, FLUXCALERR, MAG = [], [], [], []
        for i,row in enumerate(obs):
            if row[1] == band:
                MJD.append(row[0])
                FLUXCAL.append(row[3])
                FLUXCALERR.append(row[4])
                MAG.append(row[9])
        MJD, FLUXCAL, FLUXCALERR, MAG = np.array(MJD), np.array(FLUXCAL), np.array(FLUXCALERR), np.array(MAG)
        MJD = MJD - float(headerdict['SIM_PEAKMJD'])
        data = np.array([MJD, MAG]).transpose()
        data = np.array([x for x in data if x[1] != float(128)])
        if data != []:
            plt.plot(data[:,0], data[:,1], marker='x', label=band)
        # plt.plot(MJD, FLUXCAL, marker='x', label=band)
        # plt.errorbar(MJD, FLUXCAL, yerr=FLUXCALERR, capsize=3)

    plt.gca().invert_yaxis() #if plotting mag
    plt.legend()
    plt.xlabel('MJD - SIMPEAKMJD')
    plt.ylabel('MAG')
    plt.title(datfile)
    plt.show()


DIR = '../GW170817_AT2017gfo_LSST_WFD/'
for file in os.listdir(DIR):  # run through all .dat files in the specified directory
    if file.endswith(".DAT"):
        obs, headerdict = parse_observations(DIR + file)
        lightcurves(DIR + file)
