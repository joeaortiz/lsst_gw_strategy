import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


bands = ['g', 'r', 'i', 'z', 'Y']

min_ang_g = np.load('min_ang_g.npy')
min_ang_r = np.load('min_ang_r.npy')
min_ang_i = np.load('min_ang_i.npy')
min_ang_z = np.load('min_ang_z.npy')
min_ang_Y = np.load('min_ang_Y.npy')
min_ang_g = min_ang_g[~np.isnan(min_ang_g)]
min_ang_r = min_ang_r[~np.isnan(min_ang_r)]
min_ang_i = min_ang_i[~np.isnan(min_ang_i)]
min_ang_z = min_ang_z[~np.isnan(min_ang_z)]
min_ang_Y = min_ang_Y[~np.isnan(min_ang_Y)]
obs_g, obs_r, obs_i, obs_z, obs_y = 100-len(min_ang_g), 100-len(min_ang_r), 100-len(min_ang_i), 100-len(min_ang_z), 100-len(min_ang_Y)

def whisker_plot(min_ang_g, min_ang_r, min_ang_i, min_ang_z, min_ang_Y):
    data = [min_ang_g, min_ang_r, min_ang_i, min_ang_z, min_ang_Y]

    fig1 = plt.figure('boxplot_min_ang_perband')
    fig1.suptitle('Whisker plot of angular separation of nearest observation in a given band in the 24hrs following a trigger', fontsize=14, fontweight='bold')
    plt.title('100 randomly spatially and temporally located events in each band were used.'
              '\n The following number of events had no observations in the following 24hrs anywhere in the sky in the given band:'
              '\n g: %i/100, r: %i/100, i: %i/100, z: %i/100, Y: %i/100'
              '\n Horizontal line is at $sqrt(9/ \pi)$ deg as this is the angular radius of a lsst pointing (9 deg$^2$ area)' %(obs_g, obs_r, obs_i, obs_z, obs_y))
    plt.boxplot(data, labels=bands, showfliers=False)
    plt.ylabel('Angular separation (degrees)')

    for i in range(len(bands)): #add scatter of points
        y = data[i]
        x = np.random.normal(1+i, 0.04, size=len(y))
        plt.plot(x, y, 'r.', alpha=0.4)
    plt.hlines(np.sqrt(9.6/np.pi), 0, 5.5, color='green', alpha=0.6)
    plt.show()
    return


def histogram():
    fig2 = plt.figure('histogram_min_ang_perband')
    for band in bands:
        numbins = 20
        min_ang = np.load('min_ang_' + band + '.npy')
        min_ang = min_ang[~np.isnan(min_ang)]
        histdata, binedges = np.histogram(min_ang,numbins)
        binsize = binedges[1] - binedges[0]
        bincentres = binedges[:-1] + binsize/2
        #print '%i/100 of the events had no observations in the %s band in the following 24hrs anywhere in the sky' %((100-len(min_ang)), band)
        plt.plot(bincentres, histdata, label=band, alpha=0.7)

    plt.legend()
    plt.xlabel('Angular separation (degrees)')
    plt.ylabel('Count')
    fig2.suptitle('Histogram of angular separation of nearest observation in a given band in the 24hrs following a trigger', fontsize=14, fontweight='bold')
    plt.title('100 randomly spatially and temporally located events in each band were used. (nbins=%i)'
              '\n The following number of events had no observations in the following 24hrs anywhere in the sky in the given band:'
              '\n g: %i/100, r: %i/100, i: %i/100, z: %i/100, Y: %i/100' %(numbins, obs_g, obs_r, obs_i, obs_z, obs_y))
    plt.show()
    return

histogram()
whisker_plot(min_ang_g, min_ang_r, min_ang_i, min_ang_z, min_ang_Y)