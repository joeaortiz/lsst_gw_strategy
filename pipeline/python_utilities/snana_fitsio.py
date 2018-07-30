import os
import numpy as np
import pyfits
import scipy.ndimage
import pylab

pylab.ion()

# Label all of the
"""
phot = pyfits.open('DES_NONIa-01_PHOT.FITS') 
label = scipy.ndimage.label(phot[1].data['MJD'] > 0)[0] - 1
cut = (label == 2) & (phot[1].data['FLT'] == 'g ')
"""

def lightcurve(hdu, mjd_select=None):
    pylab.figure('lightcurve')
    print 'LIGHTCURVE'
    pylab.clf()
    MJD_REF = 56550
    BANDS = ['g', 'r', 'i']
    COLORS = ['blue', 'green', 'orange']
    for band, color in zip(BANDS, COLORS):
        cut_band = (hdu['FLT'].astype('a1') == band)
        pylab.plot(hdu['MJD'][cut_band] - MJD_REF, hdu['SIM_MAGOBS'][cut_band], c=color, alpha=1, lw=2)
        cut_detect = cut_band & (hdu['PHOTFLAG'] == 4096)
        cut_nondetect = cut_band & (hdu['PHOTFLAG'] != 4096)
        pylab.errorbar(hdu['MJD'][cut_band] - MJD_REF, hdu['MAG'][cut_band], yerr=hdu['MAGERR'][cut_band], ls='None', ecolor=color)
        pylab.scatter(hdu['MJD'][cut_detect] - MJD_REF, hdu['MAG'][cut_detect], c=color, edgecolor='none', label=band)
        pylab.scatter(hdu['MJD'][cut_nondetect] - MJD_REF, hdu['MAG'][cut_nondetect], c='none', edgecolor=color)
        if mjd_select:
            # SHOW THE OBSERVATION NIGHTS
            pass

    pylab.legend(loc='upper left')
    pylab.xlim(56549 - MJD_REF, 56580 - MJD_REF)
    pylab.ylim(26., 18.)
    pylab.xlabel('MJD - %i'%(MJD_REF))
    pylab.ylabel('MAG')

def selectCID(listfile, cid):
    """
    Return the photometry information for a single SN event, selected by CID
    """    
    reader = open(listfile)
    headfiles = np.char.strip(reader.readlines())
    reader.close()

    headfile_select = []
    index_select = []
    for headfile in headfiles:
        head = pyfits.open('%s/%s'%(os.path.dirname(listfile), headfile))
        index = np.nonzero(head[1].data['SNID'].astype(int) == cid)[0]
        if len(index) == 1:
            headfile_select.append(headfile)
            index_select.append(index[0])
        head.close()     
    assert len(headfile_select) == 1
    assert len(index_select) == 1

    photfile = ('%s/%s'%(os.path.dirname(listfile), headfile_select[0]))
    photfile = photfile.replace('HEAD', 'PHOT')
    phot = pyfits.open(photfile)
    label = scipy.ndimage.label(phot[1].data['MJD'] > 0)[0] - 1
    cut = (label == index_select[0])
    hdu = phot[1].data[cut]
    phot.close()
    return hdu

def search(hdu, bands=('g', 'r', 'i'), mjd_select=(56560, 56567, 56574)):
    cut = np.in1d(hdu['FLT'].astype('a1'), ('g', 'r', 'i')) & \
          np.in1d(hdu['MJD'].astype(int), mjd_select)
    n_detect_per_night = np.sum(hdu['PHOTFLAG'][cut].reshape(-1, len(bands)) == 4096, axis=1)
    detect = False
    if np.sum(n_detect_per_night >= 1) >=2:
        detect = True

    return detect
    
def fieldsView(arr, fields):
    """
    http://stackoverflow.com/questions/15182381/how-to-return-a-view-of-several-columns-in-numpy-structured-array
    """
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

#def detectChunk(headfile, bands=('g', 'r', 'i'), mjd_select=(56560, 56567, 56574)):
def detectChunk(headfile, bands=('g', 'r', 'i'), mjd_select=(56553, 56560, 56567, 56574)):
    photfile = headfile.replace('HEAD', 'PHOT')
    head = pyfits.open(headfile)
    phot = pyfits.open(photfile)

    label = scipy.ndimage.label(phot[1].data['MJD'] > 0)[0] - 1
    n_events = np.max(label) + 1

    #print label[0:1000]
    #print label[-1000:]
    #print n_events
    #print len(head[1].data)
    #print 'test', np.arange(n_events + 1)  - 0.5
    #print len(np.arange(n_events + 1)  - 0.5)

    assert n_events == len(head[1].data)

    hdu = phot[1].data
    cut = np.in1d(hdu['FLT'].astype('a1'), ('g', 'r', 'i')) & \
          np.in1d(hdu['MJD'].astype(int), mjd_select)

    label = label[cut]
    hdu = fieldsView(hdu[cut], ('MJD', 'PHOTFLAG', 'MAG', 'MAGERR'))
    
    #print np.max(label)
    #print label[0:1000]

    #return label
    #detection_criteria = (hdu['PHOTFLAG'] >= 4096) & (hdu['MAG'] < 23.) & (hdu['MAGERR'] < 0.1)
    detection_criteria = (hdu['PHOTFLAG'] >= 4096) & (hdu['MAGERR'] < 0.1)
    #n_detect_per_night = np.sum((hdu['PHOTFLAG'] >= 4096).reshape(-1, len(bands)), axis=1)
    n_detect_per_night = np.sum(detection_criteria.reshape(-1, len(bands)), axis=1)
    #label_per_night = label.reshape(len(bands), -1)[0]
    label_per_night = label.reshape(-1, len(bands)).T[0]

    #print 'max', np.max(label_per_night)
    #print n_events
    #print len(n_detect_per_night)
    #print label_per_night[0:100]
    #print n_detect_per_night[0:100]

    n_detect_per_event = np.histogram(label_per_night[n_detect_per_night >= 1], bins=np.arange(n_events + 1) - 0.5)[0] # Detected in at least 1 band

    detect_array = n_detect_per_event >= 2 # Detected on at least 2 nights
    #detect_array = detect_array & (head[1].data['SIM_PEAKMAG_r'] < 24.)
    cid_array = head[1].data['SNID'].astype(int)

    #pylab.figure()
    #pylab.scatter(head[1].data['SIM_PEAKMJD'][~detect_array], head[1].data['SIM_PEAKMAG_r'][~detect_array], c='0.5', edgecolor='none')
    #pylab.scatter(head[1].data['SIM_PEAKMJD'][detect_array], head[1].data['SIM_PEAKMAG_r'][detect_array], c='red', edgecolor='none')
    #pylab.ylim(34., 16.)
    #raw_input()

    return cid_array, detect_array

    #return label, hdu


    """
    label = scipy.ndimage.label(phot[1].data['MJD'] > 0)[0] - 1
    detect_array = np.tile(False, np.max(label))
    for ii in range(0, np.max(label)):
        print ii
        cut = (label == ii)
        hdu = phot[1].data[cut]
        detect_array[ii] = search(hdu)
    """ 
    
def detectAll(listfile):
    reader = open(listfile)
    headfiles = np.char.strip(reader.readlines())
    reader.close()
    
    cid_array = []
    detect_array = []
    for headfile in headfiles:
        cid, detect = detectChunk('%s/%s'%(os.path.dirname(listfile), headfile))
        cid_array.append(cid)
        detect_array.append(detect)

    return np.concatenate(cid_array), np.concatenate(detect_array)
