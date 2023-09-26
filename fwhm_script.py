import os

import astropy.io.fits as fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, ApertureStats
from photutils.detection import DAOStarFinder

path2data = r'D:\RoboPhot Data\кал\195'
file_list = []
dir_content = os.listdir(path2data)
fwhms = []
for ii in range(0, len(dir_content)):
    if dir_content[ii].count('.fit') or dir_content[ii].count('.fits'):
        file_list.append(path2data + '/' + dir_content[ii])
for item in file_list:
    hdu_list = fits.open(item)
    header = hdu_list[0].header
    Data = hdu_list[0].data
    hdu_list.verify('fix')
    hdu_list.close()
    mean, median, std = sigma_clipped_stats(Data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)
    sources = daofind(Data - median)
    for col in sources.colnames:
        if col not in ('id', 'npix'):
            sources[col].info.format = '%.2f'  # for consistent table output
    # sources.pprint(max_width=76)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.0)
    aperstats = ApertureStats(Data, apertures)
    fwhm = aperstats.fwhm.value * (20.1 * 60 / 1024)
    fwhms.append([np.round(min(fwhm), 2), np.round(max(fwhm), 2),
                  np.round(np.mean(fwhm), 2), np.round(np.median(fwhm), 2)])
print('min \t max \t mean \t median \t link')
for i in range(len(fwhms)):
    print(*fwhms[i], sep='\t', end='\t')
    print(file_list[i])
