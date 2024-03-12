import numpy as np
from astropy.modeling import models
import matplotlib.pyplot as plt
from scipy import ndimage
from astropy.stats import sigma_clipped_stats


##empty array
imshape = (400, 400)
image = np.zeros(imshape)

##add simple noise
Sky_noise = 20
image = image + np.random.normal(0, Sky_noise**2, imshape)
##plt.imshow(image)
##plt.show()

##add artificial stars
N_stars = 50
Amps = np.random.random(N_stars) * 50000
X = np.random.random(N_stars) * imshape[1]
Y = np.random.random(N_stars) * imshape[0]

##set target FWHM here
FWHM = 2.5
for I, A in enumerate(Amps):
    g = models.Gaussian2D(amplitude=A, x_mean=X[I], y_mean=Y[I],
                          x_stddev=FWHM/2.355, y_stddev=FWHM/2.355)
    g.render(image)
##    print(g.x_fwhm)

plt.imshow(image)
plt.show()

##apply filters
f_image = ndimage.median_filter(image, 9, mode= 'reflect')
f_image = ndimage.gaussian_filter(f_image, 3, 0, mode= 'reflect')
##calc noise etc.
mean, median, stddev = sigma_clipped_stats(f_image, sigma = 3, maxiters = 3, \
                                               cenfunc = 'median', stdfunc = 'mad_std')

##detect peaks
Peaks = f_image - (median+10*stddev)
detected_peaks = Peaks>0
labeled_im, nb_labels = ndimage.label(detected_peaks)
##check labels size
sizes = ndimage.sum(detected_peaks, labeled_im, range(nb_labels + 1))
mask_size = sizes < 5
remove_pixel = mask_size[labeled_im]
labeled_im[remove_pixel] = 0
labeled_im[labeled_im>0]=100
##redetect features
labeled_im, nb_labels = ndimage.label(labeled_im)
slices = ndimage.find_objects(labeled_im)

FWHM = []
for Slice in slices:
    ##check roundness
    X2Y = (Slice[0].stop - Slice[0].start) / (Slice[1].stop - Slice[1].start)
    if (X2Y>1.2) or (X2Y<0.8):
        continue    
    ##copy small area of the image
    Slice = image[Slice] - median
    ##index_array
    Y_index = np.arange(0, Slice.shape[0], 1)
    X_index = np.arange(0, Slice.shape[1], 1)
    ##calc centroid
    My = np.sum(Slice*Y_index[:, None])/np.sum(Slice)
    Mx = np.sum(Slice*X_index[None, :])/np.sum(Slice)

    ##calc second order moments
    Y_index = Y_index - My
    X_index = X_index - Mx
    Myy = np.sum(Slice*Y_index[:, None]*Y_index[:, None])/np.sum(Slice)
    Mxx = np.sum(Slice*X_index[None, :]*X_index[None, :])/np.sum(Slice)
    ##calc FWHM
    _fwhm = 2 * np.sqrt(0.69 * (Mxx + Myy))
##    print('Centriod: ', Mx, My, '\t FWHM: ', _fwhm)
    FWHM.append((_fwhm))
    plt.imshow(Slice)
    plt.plot(Mx, My, 'rx')
    plt.show()

print(np.median(np.asarray(FWHM)))
plt.hist(FWHM, 10)
plt.show()

