import numpy as np

def make_simlib(new_pointings_arr, simlib_file): #simlib_file is name of output file to write to

    libids = np.unique(new_pointings_arr[:,0])
    
    with open(simlib_file, 'w') as simlib:

        #write title lines
        simlib.write('SURVEY: LSST    FILTERS: ugrizY  TELESCOPE: LSST\n')
        simlib.write('USER: rbiswas     HOST: time\n')
        simlib.write('BEGIN LIBGEN\n\n')

        for libid in libids:

            pointings = new_pointings_arr[new_pointings_arr[:,0]==libid]
            n_obs = len(pointings)
            ra, decl = pointings[0,1], pointings[0,2]

            #make field header
            simlib.write('# --------------------------------------------' +'\n')
            simlib.write('LIBID:{0:11d}'.format(int(libid)) +'\n')
            simlib.write('RA: {0:+.6f}'.format(ra) + ' DECL: {0:+.6f}'.format(decl) )
            simlib.write('   NOBS: {0:11d}'.format(n_obs) + ' MWEBV: 0.01 PIXSIZE: 0.200' + '\n')
            simlib.write('#                           CCD  CCD         PSF1 PSF2 PSF2/1' +'\n')
            simlib.write('#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG' + '\n')

            #insert pointings
            #[LIBID, RA, DECL, MJD IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG]
            sorted_pointings = pointings[np.argsort(pointings[:,3])]
            bands = ['u', 'g', 'r', 'i', 'z', 'Y']

            for pointing in sorted_pointings:

                simlib.write('S: ' + '{0:5.4f} '.format(pointing[3]) + '{0:10d} '.format(int(pointing[4])))
                simlib.write('{} '.format(bands[int(pointing[5])]) + '{0:5.2f} '.format(1.))
                simlib.write("{0:5.2f} ".format(0.25) + "{0:6.2f} ".format(pointing[8]))
                simlib.write("{0:4.2f} ".format(pointing[9]) + "{0:4.2f} ".format(0.) + "{0:4.3f} ".format(0.))
                simlib.write('{0:6.2f} '.format(pointing[12]) + '{0:6.3f} '.format(pointing[13]) )
                simlib.write("{0:+7.3f}".format(-99.) + '\n')


            #make field footer
            simlib.write('END_LIBID:{0:11d}'.format(int(libid)) + '\n')
          
        simlib.write('END_OF_SIMLIB:          2293 ENTRIES')
        
    return


#new_pointings_arr = np.load('new_pointing_arrays/1_28_2_70/new_pointings_arr.npy')
#make_simlib(new_pointings_arr, '../simlibs/new.simlib')

