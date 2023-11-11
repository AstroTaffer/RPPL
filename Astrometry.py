import os

import astropy.io.fits as fits
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std
from astroquery.astrometry_net import AstrometryNet
from photutils.detection import DAOStarFinder
from astropy import units as u


def Astrometry(path, files, C):
    ast = AstrometryNet()
    ast.api_key = 'hipfhzhlzygnlvix'
    path2save = path + r'\AstrometryDone'
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
            try:
                buf = Header['CD1_1']
                hdulist.close()
                # print('done')
                continue
            except:
                pass
            Data = hdulist[0].data.copy()
            hdulist.verify('fix')
            hdulist.close()
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
                                                    radius=0.25,
                                                    downsample_factor=2,
                                                    scale_lower=1.2,
                                                    scale_upper=1.4,
                                                    scale_units='arcsecperpix'
                                                    )
            hdu = fits.PrimaryHDU(hdulist[0].data, Header + wcs_header)
            # hdulist[0].header = Header + wcs_header
            # hdulist.close()
            hdu.writeto(path2save + '\\' + name, overwrite=True)
            # hdu.writeto(item, overwrite=True)
            print('done')
        except Exception as e:
            print(e)
    return path2save


class DoAstrometry:

    def __init__(self, key):
        self.ast = AstrometryNet()
        self.ast.api_key = key

    def compute(self, input_path):
        try:
            # read file, copy data and header
            hdulist = fits.open(input_path, 'update', memmap=False)
            Header = hdulist[0].header
            try:
                buf = Header['CD1_1']
                hdulist.close()
                return True
            except:
                pass
            C = SkyCoord(f'{Header["ALPHA"]} {Header["DELTA"]}', unit=(u.hourangle, u.deg), frame='icrs')
            Data = hdulist[0].data.copy()
            hdulist.verify('fix')
            hdulist.close()
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
            # Sort sources in ascending order
            Sources.sort('flux')
            Sources.reverse()
            # ast.show_allowed_settings()
            image_width = Header['NAXIS2']
            image_height = Header['NAXIS1']
            # print(Sources)
            wcs_header = self.ast.solve_from_source_list(Sources['xcentroid'],
                                                         Sources['ycentroid'],
                                                         image_width, image_height,
                                                         solve_timeout=120,
                                                         center_ra=C.ra.degree,
                                                         center_dec=C.dec.degree,
                                                         radius=0.25,
                                                         downsample_factor=2,
                                                         scale_lower=1.2,
                                                         scale_upper=1.4,
                                                         scale_units='arcsecperpix'
                                                         )
            hdulist[0].header = Header + wcs_header
            hdulist.close()
            return True
        except Exception as e:
            print(e)
            return False
