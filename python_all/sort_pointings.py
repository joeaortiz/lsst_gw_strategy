import numpy as np

from parse_simlib import parse_simlib

libid_details, pointings = parse_simlib('../minion_1016_WFD.simlib')


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

print pointings_arr[0:10]
