import numpy as np
import matplotlib.pyplot as plt

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

#check mjd is rising
# x = np.arange(0,len(pointings_arr))
# fig3 = plt.figure()
# plt.plot(x,pointings_arr[:,3], alpha=0.5)
# plt.show()

band = -1
problems = []
mjd_dict = []
time =0
for pointing in pointings_arr:
    if time == pointing[3]:
        mjd_dict[-1][1].append(pointing[0])
        mjd_dict[-1][2].append(pointing[5])
        if pointing[5] != band:
            problems.append(len(mjd_dict))



    else:
        mjd_dict.append([pointing[3], [pointing[0]], [pointing[5]]])
    time = pointing[3]
    band = pointing[5]

print mjd_dict[0:3]
print problems
print mjd_dict[27]
print mjd_dict[244:248]


#create array with 2 columns: [overhead time, change band?] changeband = 0 for no change and 1 for change of band.
#length of overhead_times is len(mjd_band) - 1
overhead_times = np.zeros([len(pointings_arr)-1, 2])
for i in range(len(overhead_times)):
    overhead_times[i,0] = pointings_arr[i+1,0] - pointings_arr[i,0]
    time = pointings_arr[i]
    if pointings_arr[i+1,1] != pointings_arr[i,1]:
        overhead_times[i,1] = 1

print np.max(overhead_times[:,0])
print np.min(overhead_times[:,0])
print np.median(overhead_times[:,0])
print np.mean(overhead_times[:,0])
print (overhead_times[:,0] == 0).sum(), '\n'

for i in range(100):
    if overhead_times[i,0] == 0:
        print mjd_band[i+1,0], mjd_band[i,0], mjd_band[i+1,2], mjd_band[i,2]

#plot histogram
# fig2 = plt.figure()
# plt.hist(overhead_times[:,0], bins=20)
# plt.show()

def whisker_plot(data):

    fig1 = plt.figure('boxplot_min_ang_perband')
    fig1.suptitle('Whisker plot of angular separation of nearest observation in a given band in the 24hrs following a trigger', fontsize=14, fontweight='bold')
    plt.title('title')
    plt.boxplot(data, showfliers=False)
    plt.ylabel('Overhead time (days)')

    y = data
    x = np.random.normal(1, 0.04, size=len(y))
    #plt.plot(x, y, 'r.', alpha=0.1)
    plt.ylim((0,0.01))
    plt.show()
    return

# whisker_plot(overhead_times[:,0])