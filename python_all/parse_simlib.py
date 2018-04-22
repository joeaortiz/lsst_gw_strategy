import numpy as np
# from itertools import ifilter

headers = ['RA:', 'DECL:', 'NOBS:', 'MWEBV:', 'PIXSIZE:']



def parse_simlib(simlib_infile):
    #create array libid_details with each row corresponding to a different libid
    #create an array of all pointings, with [LIBID, RA, DECL] appended to the start of each pointing. Column labels after are [MJD IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG]
    with open(simlib_infile, 'r') as inp:
        libid_details= [] #list with each row giving ['LIBID:', 'RA:', 'DECL:', 'NOBS:', 'MWEBV:', 'PIXSIZE:', 'first observation index', 'last observation index']
        observations = 0
        pointings = []
        for line in inp:
            if len(line.split()) != 0:
                if line.split()[0] == 'LIBID:':
                    LIBID = line.split()[1]
                if line.split()[0] == 'RA:':
                    libid_inf = [x for x in line.split() if x not in headers]
                    libid_inf.insert(6, observations)
                    observations += int(libid_inf[2])
                    libid_inf.insert(7, observations)
                    libid_inf.insert(0, LIBID)
                    libid_details.append(libid_inf)
                if line.split()[0] == 'S:':
                    obs = line.split()[1:]
                    obs.insert(0, libid_inf[2])#insert DECL into pointings
                    obs.insert(0, libid_inf[1])#insert RA into pointings
                    obs.insert(0, libid_inf[0])#insert LIBID into pointings
                    pointings.append(obs)
                    
    bands = ['u', 'g', 'r', 'i', 'z', 'Y']

    #map bands onto numbers to we can put list of pointings into a np array and convert all strings into floats
    for row in range(len(pointings)):
        for num, band in enumerate(bands):
            if pointings[row][5] == band:
                pointings[row][5] = num
        pointings[row] = [float(i) for i in pointings[row]] #turn all entries in list into floats
        pointings[row][0], pointings[row][4], pointings[row][5] = int(pointings[row][0]), int(pointings[row][4]), int(pointings[row][5]) #some values take integers only e.g. libid


    pointings_arr = np.array(pointings)
    pointings_arr = pointings_arr[pointings_arr[:,3].argsort()] #sort in time order

    return libid_details, pointings_arr


def sayhi():
    print 'hi'
    return 

def libid_pointings(LIBID):
    #create array of each pointing for given LIBID
    pointings = []
    with open('../LSST_WFD_COADD.SIMLIB', 'r') as inp:
        copy = False
        start = 'LIBID: ' + LIBID
        end = 'END_LIBID: ' + LIBID
        for line in inp:
            if len(line.split()) != 0:
                if line.strip() == start:
                    copy = True
                elif line.strip() == end:
                    break
                elif copy:
                    pointings.append(line.split()[1:])
    return pointings[3:]

#
# libid_details, pointings = parse_simlib('../LSST_WFD_COADD.SIMLIB')
#
# print pointings[0:100]
# print len(pointings)