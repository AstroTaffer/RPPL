import os
from astropy.io import ascii
import subprocess

import astropy.io.fits as fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, ApertureStats
from photutils.detection import DAOStarFinder
import matplotlib.pyplot as plt


def sex(input_file):
    cwd = 'C:\\'
    # cwd = os.getcwd() + '\\'
    Sex = cwd + 'Sex\Extract.exe '
    dSex = ' -c ' + cwd + 'Sex\pipeline.sex'
    dPar = ' -PARAMETERS_NAME ' + cwd + 'Sex\pipeline.par'
    dFilt = ' -FILTER_NAME ' + cwd + r'Sex\tophat_2.5_3x3.conv'
    NNW = ' -STARNNW_NAME ' + cwd + 'Sex\default.nnw'

    output_file = ".".join(input_file.split('.')[:-1]) + '.cat'
    # output_file = input_file.replace('fits.gz', 'cat')

    shell = Sex + "\"" + input_file + "\"" + dSex + dPar + dFilt + NNW + ' -CATALOG_NAME ' + "\"" + output_file + "\""
    print(shell)
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    child = subprocess.run(shell, timeout=60, startupinfo=startupinfo)

    if child.returncode == 0 and os.path.isfile(output_file):
        print('Ok')
    else:
        print('Error')
        return 0, 0, 0, ''
    tbl = ascii.read(output_file)
    # indx = np.where((tbl['FWHM_IMAGE'] < 50) & (tbl['FWHM_IMAGE'] > 1))[0]
    indx = np.where((tbl['FWHM_IMAGE'] > 1) & (tbl['FLUX_ISOCOR']/tbl['FLUXERR_ISOCOR'] > 15) &
                    (tbl['FLUX_ISOCOR']/tbl['FLUXERR_ISOCOR'] < 1000) & (tbl['FLAGS'] == 0))[0]
    # & (tbl['FLAGS'] == 0) &
    # (tbl['FLUX_ISOCOR']/tbl['FLUXERR_ISOCOR'] > 15) &
    # (tbl['FLUX_ISOCOR']/tbl['FLUXERR_ISOCOR'] < 1000)
    if len(indx) < 10:
        print('Can\t find stars')
        return 0
    med_fwhm = np.round(np.median(tbl['FWHM_IMAGE'][indx]), 2)
    # med_ell = np.round(np.median(tbl['ELLIPTICITY'][indx]), 2)
    # med_bkg = np.round(np.median(tbl['BACKGROUND'][indx]), 2)
    # med_zeropoi = np.round(np.median(tbl['ZEROPOI']), 2)
    return med_fwhm


path2data = r'D:\2024\2024-01-24\Focus\FOCUS'
file_list = []
dir_content = os.listdir(path2data)
fwhms = []
for ii in range(0, len(dir_content)):
    if dir_content[ii].count('.fit') or dir_content[ii].count('.fits'):
        file_list.append(path2data + '/' + dir_content[ii])
for item in file_list:
    if '_r_' in item:
        break
    print(f'calc {item}')
    hdu_list = fits.open(item)
    header = hdu_list[0].header
    Data = hdu_list[0].data
    hdu_list.verify('fix')
    hdu_list.close()
    fwhm = sex(item)
    if fwhm == 0:
        continue
    # mean, median, std = sigma_clipped_stats(Data, sigma=3.0)
    # daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)
    # sources = daofind(Data - median)
    # try:
    #     for col in sources.colnames:
    #         if col not in ('id', 'npix'):
    #             sources[col].info.format = '%.2f'  # for consistent table output
    # except:
    #     continue
    # # sources.pprint(max_width=76)
    # positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    # apertures = CircularAperture(positions, r=4.0)
    # aperstats = ApertureStats(Data, apertures)
    # fwhm = aperstats.fwhm.value * (22 * 60 / 1024)
    # fwhms.append([np.round(min(fwhm), 2), np.round(max(fwhm), 2),
    #               np.round(np.mean(fwhm), 2), np.round(np.median(fwhm), 2)])
    fwhms.append([header['FOCUS'], np.round(np.median(fwhm), 2)])

fwhms = np.array(fwhms)
plt.plot(fwhms[:, 0], fwhms[:, 1], 'b.')
plt.savefig('focus.pdf')
np.savetxt('focus.txt', fwhms)
# print('min \t max \t mean \t median \t link')
# for i in range(len(fwhms)):
#     print(*fwhms[i], sep='\t', end='\t')
#     print(file_list[i])


