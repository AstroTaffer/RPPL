import os

import astropy.io.fits as fits
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.stats import mad_std
from astroquery.astrometry_net import AstrometryNet
from photutils.detection import DAOStarFinder


def Astrometry(path, files, C):
    ast = AstrometryNet()
    ast.api_key = 'hipfhzhlzygnlvix'
    path2save = path + r'\done'
    if not os.path.exists(path2save):
        os.mkdir(path2save)
    for count, item in enumerate(files):
        try:
            name = item.split('\\')[-1]
            print(f"Astrometry ({count}/{len(files)}): working on {name}")
            # read file, copy data and header
            hdulist = fits.open(item, 'update', memmap=False)
            # del hdulist[0].header['COMMENT']
            Header = hdulist[0].header
            Data = hdulist[0].data.copy()
            # hdulist.verify('fix')
            # gaussian convolution
            kernel = Gaussian2DKernel(x_stddev=1)
            Data = convolve(Data, kernel)
            # extract background
            Data -= np.median(Data)
            Bkg_sigma = mad_std(Data)
            # # mask bad row
            # mask = np.zeros(Data.shape, dtype=bool)
            # mask[90:110, 0:Data.shape[1]] = True
            daofind = DAOStarFinder(fwhm=4.5, threshold=5. * Bkg_sigma, sharplo=0.25)
            # Sources = daofind(Data, mask=mask)
            Sources = daofind(Data)
            # print(Sources.info)
            # plt.imshow(Data, cmap=cm.Greys_r, aspect='equal',
            #            norm=Normalize(vmin=-30, vmax=150), interpolation='nearest')
            # plt.scatter(Sources['xcentroid'], Sources['ycentroid'], s=40, facecolors='none', edgecolors='r')
            # plt.show()
            # Sort sources in ascending order
            Sources.sort('flux')
            Sources.reverse()
            # ast.show_allowed_settings()
            image_width = Header['NAXIS2']
            image_height = Header['NAXIS1']
            # print(Sources)
            wcs_header = ast.solve_from_source_list(Sources['xcentroid'],
                                                    Sources['ycentroid'],
                                                    image_width, image_height,
                                                    solve_timeout=120,
                                                    center_ra=C.ra.degree,
                                                    center_dec=C.dec.degree,
                                                    radius=0.1,
                                                    downsample_factor=2,
                                                    scale_lower=0.6,
                                                    scale_upper=1.5,
                                                    scale_units='arcsecperpix'
                                                    )
            # hdu = fits.PrimaryHDU(hdulist[0].data, Header + wcs_header)
            hdulist[0].header = Header + wcs_header
            hdulist.close()
            # hdu.writeto(path2save + '\\' + name, overwrite=True)
            # hdu.writeto(item, overwrite=True)
            print('done')
        except Exception as e:
            print(e)
    return path

# paths = [r'D:\RoboPhot Data\Raw Images\2023_09 FLI\2023_09_06 GSC2314–0530\i',
#          r'D:\RoboPhot Data\Raw Images\2023_09 FLI\2023_09_06 GSC2314–0530\r']
# coords = '02:20:50.9 +33:20:46.6'
# C = SkyCoord(coords, unit=(u.hourangle, u.deg), frame='icrs')
# Astrometry(r'C:\Users\RoboPhot_win\Desktop\test', C)
# for path in paths:
#     C = SkyCoord(coords, unit=(u.hourangle, u.deg), frame='icrs')
#     Astrometry(path, C)
