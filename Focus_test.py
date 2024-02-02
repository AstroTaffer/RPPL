import os

import numpy as np
from astropy.io import fits as pyfits
from astropy.stats import sigma_clipped_stats  # , sigma_clip
from matplotlib import pyplot as plt
from scipy import ndimage
from glob import glob

Path = r'D:\2024\2024-01-24\Focus\FOCUS\\*.fits'
List = glob(Path)

Focus = []
for file in List:

    hdulist = pyfits.open(file)
    Header = hdulist[0].header
    Data_0 = hdulist[0].data
    hdulist.close()
    print(file, Header['FOCUS'])

    ##    mean, median, stddev = sigma_clipped_stats(Data_0, sigma = 3, maxiters = 3, \
    ##                                               cenfunc = 'median', stdfunc = 'mad_std')

    Data = ndimage.median_filter(Data_0, 9, mode='reflect')
    Data = ndimage.gaussian_filter(Data, 3, 0, mode='reflect')

    mean, median, stddev = sigma_clipped_stats(Data, sigma=3, maxiters=3, \
                                               cenfunc='median', stdfunc='mad_std')

    ##search features
    Peaks = Data - (median + 10 * stddev)
    detected_peaks = Peaks > 0
    labeled_im, nb_labels = ndimage.label(detected_peaks)

    sizes = ndimage.sum(detected_peaks, labeled_im, range(nb_labels + 1))
    mask_size = sizes < 10
    remove_pixel = mask_size[labeled_im]
    labeled_im[remove_pixel] = 0
    labeled_im[labeled_im > 0] = 100
    labeled_im, nb_labels = ndimage.label(labeled_im)
    slices = ndimage.find_objects(labeled_im)
    FWHM = []
    for Slice in slices:
        X2Y = (Slice[0].stop - Slice[0].start) / (Slice[1].stop - Slice[1].start)
        if (X2Y > 1.2) or (X2Y < 0.8):
            continue
        Slice = Data_0[Slice] - median
        X_index = np.arange(0, Slice.shape[1], 1)
        Y_index = np.arange(0, Slice.shape[0], 1)
        Mxx = np.sum(Slice * X_index[None, :] * X_index[None, :]) / np.sum(Slice)
        Myy = np.sum(Slice * Y_index[:, None] * Y_index[:, None]) / np.sum(Slice)
        FWHM.append(2 * np.sqrt(0.69 * (Mxx + Myy)))
    ##        print(2 * np.sqrt (0.69 * (Mxx + Myy)))
    ##        plt.imshow(Slice)
    ##        plt.show()
    Focus.append((np.median(FWHM), Header['FOCUS']))
##    plt.hist(FWHM, 10)
##    plt.show()
##    plt.imshow(labeled_im)
##    plt.show()

Focus = np.asarray(Focus)
plt.plot(Focus[:, 1], Focus[:, 0], 'b.')
plt.show()

##https://astro.uni-bonn.de/~sysstw/lfa_html/iraf/images.tv.imexamine.html
##Mxx = sum (x * x * I) / sum (I)
##Myy = sum (y * y * I) / sum (I)
##Mxy = sum (x * y * I) / sum (I)
##r = 2 * sqrt (ln (2) * (Mxx + Myy))
##e = sqrt ((Mxx - Myy) ** 2 + (2 * Myy) ** 2)
