import json
from os.path import exists
import os
from sqlalchemy import create_engine
from sqlalchemy import URL
import astropy.io.fits as fits
import numpy as np
from photutils.background import MedianBackground, Background2D
from scipy import ndimage
from astropy.stats import sigma_clipped_stats, SigmaClip


def read_fits_file(file_name):
    with fits.open(file_name) as hdu_list:
        hdu_list.verify("fix")
        header = hdu_list[0].header
        data = hdu_list[0].data
    return header, data


def check_out_directory(out_dir):
    if not exists(out_dir):
        os.makedirs(out_dir)


# LAST UPDATE 05-09-2023
def restore_default_config():
    settings = {
        "IN_DIRS": "",
        "CALIBRATION_TYPE": "DO_BOTH",
        "OUT_DIR": "",
        "EPSILON_TEMP": 1,
        "MIN_IMAGES_FOR_MASTER": 5,
        "CALIBRATED_BITPIX": 16}
    with open(f"default_config.json", "w") as confile:
        json.dump(settings, confile, indent=4)


def get_fwhm_data(input_file):
    try:
        with fits.open(input_file, memmap=False) as hdulist:
            h = hdulist[0].header
            fwhm = float(h['FWHM'])
            ell = float(h['ELL'])
            # nstars = h['NSTARS']
            bkg = float(h['BKG'])
            if 'nan' in fwhm:
                return 0, 0, 0
        return fwhm, ell, bkg
    except:
        try:
            return calc_fwhm(input_file)
        except:
            return 0, 0, 0


def connect_to_db():
    # DB params
    db_host = "192.168.240.5"
    db_name = "postgres"
    db_user = "postgres"
    db_pass = 'Vjlekm_hfccnjzybz'
    url_object = URL.create(
        "postgresql",
        username=db_user,
        password=db_pass,
        host=db_host,
        database=db_name,
    )
    try:
        # open connection
        return create_engine(url_object)
    except Exception as e:
        print('Error while connecting to db')
        print(e)


def make_out_path(out_frame_fp):
    path2save = "\\".join(out_frame_fp.split('\\')[:-1])
    if not os.path.exists(path2save):
        os.makedirs(path2save)


def calc_fwhm(path):
    with fits.open(path, memmap=False, mode='update') as hdulist:
        header = hdulist[0].header
        image = hdulist[0].data.copy()
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(image, (50, 50), filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        # apply filters
        f_image = ndimage.median_filter(image, 9, mode='reflect')
        f_image = ndimage.gaussian_filter(f_image, 3, 0, mode='reflect')
        # calc noise etc.
        mean, median, stddev = sigma_clipped_stats(f_image, sigma=3, maxiters=3,
                                                   cenfunc='median', stdfunc='mad_std')
        # detect peaks
        Peaks = f_image - (median + 10 * stddev)
        detected_peaks = Peaks > 0
        labeled_im, nb_labels = ndimage.label(detected_peaks)
        # check labels size
        sizes = ndimage.sum(detected_peaks, labeled_im, range(nb_labels + 1))
        mask_size = sizes < 5
        remove_pixel = mask_size[labeled_im]
        labeled_im[remove_pixel] = 0
        labeled_im[labeled_im > 0] = 100
        # redetect features
        labeled_im, nb_labels = ndimage.label(labeled_im)
        slices = ndimage.find_objects(labeled_im)

        FWHM = []
        ELL = []
        for Slice in slices:
            # check roundness
            X2Y = (Slice[0].stop - Slice[0].start) / (Slice[1].stop - Slice[1].start)
            if (X2Y > 1.2) or (X2Y < 0.8):
                continue
                # copy small area of the image
            Slice = image[Slice] - median
            # index_array
            Y_index = np.arange(0, Slice.shape[0], 1)
            X_index = np.arange(0, Slice.shape[1], 1)
            # calc centroid
            My = np.sum(Slice * Y_index[:, None]) / np.sum(Slice)
            Mx = np.sum(Slice * X_index[None, :]) / np.sum(Slice)

            # calc second order moments
            Y_index = Y_index - My
            X_index = X_index - Mx
            Myy = np.sum(Slice * Y_index[:, None] * Y_index[:, None]) / np.sum(Slice)
            Mxx = np.sum(Slice * X_index[None, :] * X_index[None, :]) / np.sum(Slice)
            # calc FWHM
            M = Mxx + Myy
            _fwhm = 2 * np.sqrt(0.69 * M)

            # sigmax = np.sqrt(Mxx)
            # sigmay = np.sqrt(Myy)
            ell = 1 - min([Mxx, Myy]) / M

            #     print('Centriod: ', Mx, My, '\t FWHM: ', _fwhm)
            FWHM.append(_fwhm)
            ELL.append(ell)
        fwhm = np.round(np.nanmedian(np.asarray(FWHM)), 1)
        ell = np.round(np.nanmedian(np.asarray(ELL)), 2)
        stars_num = len(FWHM)
        b = np.round(bkg.background_median, 2)
        fwhm_card = fits.Card('FWHM', 'nan' if np.isnan(fwhm) else fwhm, 'Median FWHM')
        ell_card = fits.Card('ELL', 'nan' if np.isnan(ell) else ell, 'Median ellipticity')
        stars_card = fits.Card('NSTARS', 'nan' if np.isnan(stars_num) else stars_num, "Stars on frame")
        bkg_card = fits.Card('BKG', 'nan' if np.isnan(b) else b, "Median background")
        hdulist[0].header.append(fwhm_card)
        hdulist[0].header.append(ell_card)
        hdulist[0].header.append(stars_card)
        hdulist[0].header.append(bkg_card)
    if np.isnan(fwhm):
        return 0, 0, 0
    return fwhm, ell, b
