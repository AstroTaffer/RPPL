import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from matplotlib import pyplot as plt
from numpy import arange
from photutils.aperture import CircularAperture
from photutils.centroids import centroid_sources, centroid_com
from photutils.utils import CutoutImage


def get_com(Data, Cat, Bbox):
    x, y = centroid_sources(Data, Cat['X'], Cat['Y'], box_size=Bbox,
                            centroid_func=centroid_com)
    Cat['X'] = x
    Cat['Y'] = y
    new_Max = []
    for Obj in Cat:
        new_Max.append(CutoutImage(Data, (Obj['Y'], Obj['X']),
                                   (5, 5), mode='partial').data.max())
    Cat.add_column(new_Max, name='Max')
    return Cat


def Get_GAIA(RA, DEC, R, V_lim, Cat_len):
    My_Cat = Table()
    c = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree, frame='icrs')
    V = Vizier(columns=['RAJ2000', 'DEJ2000',
                        'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag',
                        'RPmag', 'e_RPmag', "+_r"],
               column_filters={'BPmag': '>0', 'BPmag': '<' + str(V_lim)})
    V.ROW_LIMIT = Cat_len
    Vizier_result = V.query_region(c, radius=Angle(R * u.arcmin), catalog=['I/350/gaiaedr3'])  # 'I/345'
    if len(Vizier_result) != 0:
        Vizier_stars = Vizier_result[0]
        #         print(Vizier_stars.info())
        My_Cat['ID'] = arange(0, len(Vizier_stars), 1, 'int16')
        My_Cat['Ra'] = Vizier_stars['RAJ2000']
        My_Cat['Dec'] = Vizier_stars['DEJ2000']
        My_Cat['Dist'] = Vizier_stars['_r']
        My_Cat['B'] = Vizier_stars['BPmag']
        My_Cat['R'] = Vizier_stars['RPmag']
        My_Cat['G'] = Vizier_stars['Gmag']
    else:
        print('Vizier result is empty')

    #     ascii.write(My_Cat, 'My_Cat.txt', overwrite=True, delimiter='\t', format='commented_header')
    return My_Cat


def GetTrend(Flux, Cat):
    cat = Cat.copy()

    # delete bad stars
    Index = np.where(np.isnan(Flux)[0])
    if 0 in Index:
        print('Object has NaNs!')
        return 1, 0, 0
    Flux = np.delete(Flux, Index, axis=1)
    cat.remove_rows(Index)

    # main circle
    while True:
        Ensemble = Flux[:, 1:]
        Trend = np.sum(Ensemble, axis=1)
        Trend = Trend / np.mean(Trend)
        Ensemble = Ensemble / Trend[:, np.newaxis]

        # find and remove the worst star
        Std = np.std(Ensemble, 0) / np.sqrt(np.mean(Ensemble))
        if np.max(Std) > 3:
            Index = np.argmax(Std) + 1
            print('Delete star #', cat['ID'][Index],
                  ' with STD=', np.max(Std))

            Flux = np.delete(Flux, Index, axis=1)
            cat.remove_rows(Index)
        else:
            print('Stars in ensemble: ', Flux.shape[1])
            break
    return Trend, cat


def DrawMap(Image, Size, Header, Cat, Name, RAper):
    wcs = WCS(Header)
    Image = np.log10(Image)
    X = Image.shape[1] / 2
    Y = Image.shape[0] / 2
    _mean, _median, _std = sigma_clipped_stats(Image[Y - 50:Y + 50, X - 50:X + 50])
    _max = _median + 10. * _std
    _min = _median - 1. * _std

    plt.switch_backend('pdf')  # для фикса какой-то тупой ошибки в нарнии pyplot
    fig = plt.figure(figsize=(7, 7))

    ax = plt.subplot(projection=wcs, position=[0.1, 0.1, 0.8, 0.8])
    plt.imshow(Image, vmin=_min, vmax=_max, cmap='gray_r')

    ax.set_xlim(int((Header['Naxis2'] - Size) / 2),
                int((Header['Naxis2'] + Size) / 2))
    ax.set_ylim(int((Header['Naxis1'] - Size) / 2),
                int((Header['Naxis1'] + Size) / 2))

    # XY = wcs.all_world2pix(Cat['Ra'], Cat['Dec'], 0)
    # XY = np.vstack((XY[0], XY[1])).T
    XY = np.vstack((Cat["X"], Cat["Y"])).T
    aper = CircularAperture(XY, r=RAper)
    aper.plot(color='blue', lw=1.5, alpha=0.5)
    aper = CircularAperture(XY[0], r=RAper)
    aper.plot(color='red', lw=1.5, alpha=0.8)
    for i, txt in enumerate(Cat['ID']):
        plt.annotate(txt, (XY[i, 0], XY[i, 1]), color='blue', alpha=0.8)

    if Header['CD2_2'] < 0:
        plt.gca().invert_xaxis()
    if Header['CD1_1'] > 0:
        plt.gca().invert_yaxis()

    # Title = Object + ', ' + DT + '\n'
    try:
        Title = f'Filter={Header["FILTER"]}'
    except:
        Title = 'Filter=V'
    Title += ', aperture radius =' + '{:.1f}'.format(3600 * RAper * np.sqrt(Header['CD1_1']**2+Header['CD1_2']**2))+'"'
    plt.title(Title)
    ax.coords[1].set_ticklabel(rotation=90)
    ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm:ss')
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')
    ax.coords.grid(color='blue', ls='--', alpha=0.7)
    # plt.show()
    fig.savefig(Name)

# FastPlot(r'C:\Users\User\Desktop\2023_08_21 GSC2314–0530')
