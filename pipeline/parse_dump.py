"""Create an array of the simulated KN"""

import os
import numpy as np


def findall_KN(): #return [MJD, RA, DECL]

    #move to directory containing dump file
    os.chdir('/data/des41.b/data/SNDATA_ROOT/SIM/GW170817_AT2017gfo_LSST_WFD/')

    #count number of KN simulations
    dump = open('GW170817_AT2017gfo_LSST_WFD.DUMP')
    numsims = len(dump.readlines()) - 10
    print 'Number of simulations: ', numsims

    #fill all_KN matrix
    all_KN = np.zeros([numsims, 3])
    i=0
    with open('GW170817_AT2017gfo_LSST_WFD.DUMP', 'r') as inp:
        for line in inp:
            if len(line.split()) != 0:
                if line.split()[0] == 'SN:':
                    all_KN[i,:] = [float(line.split()[10]), float(line.split()[6]), float(line.split()[7])]
                    i+=1
            
    return all_KN
