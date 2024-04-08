import os
import astropy.io.fits as fits
import astropy.wcs as wcs
from astropy.stats import SigmaClip
from astropy.time import Time
from donuts import Donuts
from photutils.background import Background2D, MedianBackground
from photutils.aperture import aperture_photometry
from astropy.io import ascii

from Astrometry import Astrometry
from Photometry_Utils import *

# '02:20:50.9 +33:20:46.6' GSC2314–0530
# RA/Dec (J2000, epoch 2015.5): 20:00:16.06 +22:09:43.54 TIC 424112727
# RA/Dec (J2000, epoch 2015.5): 21:07:55.93 +48:34:56.27 TIC 357872559


def Core(path2data, do_wcs=False):
    Cat_R = 10  # catalog each cone radius
    Cat_len = 400  # catalog max length
    V_lim = 15  # catalog max depth
    Bbox = 9  # centering/morphology box half-size (pixels)
    Raper = 10  

    # read directory and create list of fits-files
    file_list = []
    for f in os.listdir(path2data):
        if f.count('.fits') or f.count('.fit'):
            file_list.append(path2data + '\\' + f)

    hdu = fits.open(file_list[0])
    header = hdu[0].header
    date = header['DATE-OBS'].split('T')[0]
    filter = header['FILTER']
    obj = header['OBJNAME']
    # obj = 'GSC2314-0530'
    coords = f'{header["ALPHA"]} {header["DELTA"]}'
    # coords = '02:20:50.9 +33:20:46.6'
    hdu.close()
    C = SkyCoord(coords, unit=(u.hourangle, u.deg), frame='icrs')  # , obstime='J2015.5'

    # do astrometry
    if do_wcs:
        path2data = Astrometry(path2data, file_list, C)

    file_list = []
    for f in os.listdir(path2data):
        if f.count('.fits') or f.count('.fit'):
            file_list.append(path2data + '\\' + f)

    # set paths and suffix for output files
    Path2Save = path2data + '/Photometry'
    if not os.path.exists(Path2Save):
        os.makedirs(Path2Save)

    # clean old photometry
    for f in os.listdir(Path2Save):
        if f.count('Phot'):
            os.remove(Path2Save+'/'+f)

    print('Download GAIA')
    Catalog = Get_GAIA(C.ra.degree, C.dec.degree, Cat_R, V_lim, Cat_len)
    ascii.write(Catalog, Path2Save + '/Cat.txt', overwrite=True, delimiter='\t')

    # open log file
    df = open(Path2Save + '/Time.txt', 'w')
    # df.write('JD\tEXPTIME\t' +
    #          'Sky\tX\tY\tMax\tShiftXDon\tShiftYDon\tShiftXCom\tShiftYCom\n')
    df.write('JD\t' +
             'EXPTIME\tEXTINCT\t' +
             'SKY-TEMP\t' +
             'Sky\tX\tY\tMax\tShiftXDon\tShiftYDon\tShiftXCom\tShiftYCom\n')
    output = open(Path2Save + '/Phot' + '.txt', 'a')
    counter = 0
    shifts_don = []
    shifts_com = []
    don = None
    obj_xy = None
    for f in file_list:
        counter = counter + 1
        print(f'----------<Frame: {counter}/{len(file_list)}>-----------')
        print(f.split('/')[-1])
        hduList = fits.open(f)
        Header = hduList[0].header
        Data = hduList[0].data
        hduList.verify('fix')
        hduList.close()
        w = wcs.WCS(Header)

        jd = Time(Header['DATE-OBS'], scale='utc').jd

        X, Y = w.all_world2pix(Catalog['Ra'], Catalog['Dec'], 0)
        Catalog.add_columns([X, Y], names=['X', 'Y'])
        indx = np.where((Catalog['X'] < 2048/Header['XBINNING']) & (Catalog['Y'] < 2048/Header['YBINNING']))[0]
        Catalog = Catalog[indx]

        # delete background
        sigmaclip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(Data, (100, 100), filter_size=(9, 9),
                           sigma_clip=sigmaclip, bkg_estimator=bkg_estimator)
        Data_without_background = Data - bkg.background
        Sky = np.median(bkg.background)
        print('Sky={0:.1f}'.format(Sky))

        # centroiding and image properties
        try:
            Catalog_COM = get_com(Data_without_background, Catalog, Bbox)  # center of mass / fast and better!
            #     Catalog = get_1dg(Data, Catalog, Bbox) # Gauss 1d fit
        except Exception as ex:
            print(ex)
            Catalog.remove_columns(names=('X', 'Y'))
            continue

        if counter == 1:
            DrawMap(Data, 256, Header, Catalog_COM, Path2Save + r'\field_com.pdf', Raper)
            don = Donuts(refimage=f, image_ext=0, overscan_width=24, prescan_width=24,
                         border=64, normalise=True, exposure='EXPTIME', subtract_bkg=True, ntiles=32)
            shifts_don = [0, 0]
            shifts_com = [0, 0]
            obj_xy = [Catalog_COM['X'][0], Catalog_COM['Y'][0]]
        else:
            shift = don.measure_shift(f)
            shifts_don = [shift.x.value, shift.y.value]
            shifts_com = [obj_xy[0] - Catalog_COM['X'][0], obj_xy[1] - Catalog_COM['Y'][0]]

        df.write('{:.7f}'.format(jd) + '\t')  # JD
        # df.write('{:.7f}'.format(Header['jd']) + '\t')  # JD
        df.write('{:.1f}'.format(Header['EXPTIME']) + '\t')  # EXPTIME

        df.write('{:.2f}'.format(Header['EXTINCT']) + '\t' if Header['EXTINCT'] != 'UNKNOWN' else '\t')
        df.write('{:.2f}'.format(Header['SKY-TEMP']) + '\t' if Header['SKY-TEMP'] != 'UNKNOWN' else '\t')
        df.write('{:.2f}'.format(Sky) + '\t')  # Sky
        df.write('{0:.3f}\t{1:.3f}\t{2:.1f}\t'.format(Catalog_COM['X'][0],
                                                      Catalog_COM['Y'][0],
                                                      Catalog_COM['Max'][0]))
        df.write(str(np.round(shifts_don[0], 2)) + '\t')
        df.write(str(np.round(shifts_don[1], 2)) + '\t')
        df.write(str(np.round(shifts_com[0], 2)) + '\t')
        df.write(str(np.round(shifts_com[1], 2)) + '\n')
        print('Max count={0:.1f}'.format(Catalog_COM['Max'][0]))

        Positions = np.transpose([Catalog_COM['X'], Catalog_COM['Y']])
        Stellar_aper = CircularAperture(Positions, r=Raper)
        Stellar_phot = aperture_photometry(Data_without_background, Stellar_aper, method='exact')

        # write to log
        Flux = Stellar_phot['aperture_sum'].value
        np.savetxt(output, Flux, fmt='%1.3f', delimiter='\t', newline='\t')
        output.write('\n')

        Catalog_COM.remove_columns(names=('X', 'Y', 'Max'))
    output.close()
    df.close()
    Plot_Curve(Path2Save, obj, date, filter)


def Plot_Curve(path2phot, objname, date, filter):
    # objname = 'GSC2314-0530'
    # date = '2023_09_06'
    camera = f'FLI filter {filter}'
    # path2phot = path2data + r'\Photometry'
    flux = np.genfromtxt(path2phot + r'\Phot.txt')
    cat = ascii.read(path2phot + r'\Cat.txt')
    time = ascii.read(path2phot + r'\Time.txt')
    Trend, Cat = GetTrend(flux, cat)
    Flux_Corr = flux / Trend[:, np.newaxis]

    # Condition_Report(objname, time, filter, path2phot)

    # tar_Flux = Flux_Corr[:, 0]
    m = -2.5 * np.log10(Flux_Corr) + 2.5 * np.log10(time['EXPTIME'])
    s = np.nanstd(m, axis=0)
    err = np.round(np.sqrt(flux[:, 0] + 3.14 * 5 ** 2 * (time['Sky'])) / (flux[:, 0]), 5)
    ZERO = int(np.max(time['JD']))
    t = time['JD'] - ZERO

    fig, axs = plt.subplots(3, 1, figsize=(6, 7), dpi=125)
    fig.suptitle(f'{objname}, {date}, camera: {camera}', fontsize=8)
    # flux vs JD
    left = 0.11
    pos = [left, 0.71, 0.85, 0.2]
    axs[0].set_position(pos)
    axs[0].errorbar(t, m[:, 0], err, fmt='r.', label=r'1$\sigma$ errorbar, $\bar\sigma=$' + str(np.round(np.nanmean(err), 5)), markersize=3)
    axs[0].legend(loc=0, fontsize=6)
    axs[0].set_ylabel('Instrumental mag', fontsize=6)
    axs[0].set_xlabel(f'JD - {ZERO}', fontsize=6)
    axs[0].set_title(f'{objname}', fontsize=6)
    locs = axs[0].get_xticks()
    tic = Time(locs, format='jd')
    x_ticks_labels = []
    for x in tic:
        # x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
        x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
    axs[0].set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    axs[0].tick_params(axis='both', labelsize=6, direction='in')
    axs[0].invert_yaxis()
    axs[0].grid()

    pos = [left, 0.35, 0.54, 0.28]
    axs[1].set_position(pos)
    axs[1].plot(t, Flux_Corr[:, 0] / np.median(Flux_Corr[:, 0]), 'r*',
                label=f'{objname}, BP=' + '{:.2f}'.format(cat['B'][0]) +
                      r', $\sigma$=' + '{:.4f}'.format(np.std(Flux_Corr[:, 0] / np.median(Flux_Corr[:, 0]))))
    # plot 6 best stars
    S_Flux = Flux_Corr[:, 1:] / np.nanmedian(Flux_Corr[:, 1:], axis=0)
    IDs = cat['ID'][1:]
    STD = np.nanstd(S_Flux, 0)

    Sort = np.argsort(STD)
    S_Flux = S_Flux[:, Sort]
    IDs = IDs[Sort]
    # plot reference stars

    for i, Val in enumerate(IDs[0:7]):
        axs[1].plot(t, S_Flux[:, i] - 0.05 * (i + 1), '.',
                    label='Ref#' + str(Val) + ', BP=' + '{:.2f}'.format(cat['B'][Val]) +
                          r', $\sigma$=' + '{:.4f}'.format(np.std(S_Flux[:, i])))
    axs[1].legend(fontsize=6, bbox_to_anchor=(1.005, 1.))
    axs[1].set_ylabel('Normalized flux', fontsize=5)
    axs[1].set_xlabel(f'JD-{ZERO}', fontsize=6)
    axs[1].tick_params(axis='both', labelsize=6, direction='in')
    axs[1].set_title('Reference stars', loc='left', fontsize=6)
    axs[1].grid()

    pos = [left, 0.07, 0.85, 0.2]
    axs[2].set_position(pos)
    axs[2].plot(time['JD']-ZERO, time['EXTINCT'], 'b.')
    axs[2].set_ylabel('extinction (mag)', fontsize=6)
    axs[2].set_xlabel(f'JD-{ZERO}', fontsize=6)
    axs[2].tick_params(axis='both', labelsize=6, direction='in')
    # axs[2].set_title('std vs magnitudes', loc='left', fontsize=6)
    # axs[2].set_ylim(0, 0.05)
    axs[2].grid()
    fig.savefig(rf'{path2phot}\plot_{objname}_{date}_{camera}.pdf')

    shift_fig, shift_ax = plt.subplots(3, 1, figsize=(6, 7), dpi=125)
    left = 0.1
    width = 0.85
    high = 0.18
    pos = [left, 0.7, width, high]
    shift_ax[0].set_position(pos)
    shift_ax[0].plot(t, time['ShiftXDon'], 'r.', label='X shift')
    shift_ax[0].plot(t, time['ShiftYDon'], 'g.', label='Y shift')
    shift_ax[0].set_title('Donuts shifts')
    shift_ax[0].grid()
    shift_ax[0].set_xlabel(f'JD-{ZERO}')
    shift_ax[0].set_ylabel('pixel')
    shift_ax[0].legend()

    pos = [left, 0.4, width, high]
    shift_ax[1].set_position(pos)
    shift_ax[1].plot(t, time['ShiftXCom'], 'r.', label='X shift')
    shift_ax[1].plot(t, time['ShiftYCom'], 'g.', label='Y shift')
    shift_ax[1].set_title('Com shifts')
    shift_ax[1].grid()
    shift_ax[1].set_xlabel(f'JD-{ZERO}')
    shift_ax[1].set_ylabel('pixel')
    shift_ax[1].legend()

    pos = [left, 0.1, width, high]
    shift_ax[2].set_position(pos)
    shift_ax[2].plot(t, time['ShiftXCom'] - time['ShiftXDon'], 'r.', label='X')
    shift_ax[2].plot(t, time['ShiftYCom'] - time['ShiftYDon'], 'g.', label='Y')
    shift_ax[2].set_title('Delta shifts')
    shift_ax[2].grid()
    shift_ax[2].set_xlabel(f'JD-{ZERO}')
    shift_ax[2].set_ylabel('pixel')
    shift_fig.savefig(rf'{path2phot}\plot_shifts.pdf')


Core(r'E:\GPX-TF16E-48\robophot\2024-03-06\i', False)
# Core(r'D:\2024\2024-01-22\GPX-TF16E-48\CALIBRATED\r', True)
# Core(r'D:\2023\2023-12-13\GPX-TF16E-48\CALIBRATED\r', False)
# Plot_Curve(r'C:\Users\User\Desktop\Tempo\2023_09_06 GSC2314-0530\DO_BOTH\i\Photometry',
#            'GSC2314–0530', '2023-09-06', 'i')
# Plot_Curve(r'C:\Users\User\Desktop\Tempo\2023_09_06 GSC2314-0530\DO_BOTH\r\Photometry',
#            'GSC2314–0530', '2023-09-06', 'r')
