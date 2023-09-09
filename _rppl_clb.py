from warnings import warn

import numpy as np
import astropy.io.fits as fits
from astropy.table import unique, Table, vstack
from astropy.stats import sigma_clip

from utils import read_fits_file, check_out_directory


class _Calibrator:
    def __init__(self):
        self.settings = None
        self.images_table = None

    def make_master_dark(self):
        check_out_directory(self.settings["OUT_DIR"])

        imgs_type = []
        imgs_binning = []
        imgs_filter = []
        imgs_exptime = []
        imgs_ccdtemp = []
        imgs_filepath = []

        sel_table_0 = self.images_table[np.where(((self.images_table["IMAGETYP"] == "Dark") |
                                                  (self.images_table["IMAGETYP"] == "Dark Frame")))]
        if len(sel_table_0) == 0:
            print("Can't create master dark frames - No dark images found\n")
            return

        bin_factors = unique(sel_table_0, "BINNING")["BINNING"]
        for _ in bin_factors:
            sel_table_1 = sel_table_0[np.where(sel_table_0["BINNING"] == _)]
            filters = unique(sel_table_1, "FILTER")["FILTER"]
            for __ in filters:
                sel_table_2 = sel_table_1[np.where(sel_table_1["FILTER"] == __)]
                exptimes = unique(sel_table_2, "EXPTIME")["EXPTIME"]
                for ___ in exptimes:
                    sel_table_3 = sel_table_2[np.where(sel_table_2["EXPTIME"] == ___)]
                    mean_ccd_temp = np.round(np.mean(sel_table_3["CCD-TEMP"]), 3)
                    sel_table_4 = sel_table_3[np.where(abs(sel_table_3["CCD-TEMP"] - mean_ccd_temp) <=
                                                       self.settings["EPSILON_TEMP"])]
                    if len(sel_table_4) < self.settings["MIN_IMAGES_FOR_MASTER"]:
                        warn(f"Less than {self.settings['MIN_IMAGES_FOR_MASTER']} dark images with BINNING {_}, "
                             f"FILTER {__}, EXPTIME {___} and CCD-TEMP near {mean_ccd_temp} - skipped")
                        continue
                    print(f"Found {len(sel_table_4)} dark images with BINNING {_}, "
                          f"FILTER {__}, EXPTIME {___} and CCD-TEMP near {mean_ccd_temp}")

                    buff_hdr = read_fits_file(sel_table_4["FILEPATH"][0])[0]
                    dark_images = np.zeros((len(sel_table_4), buff_hdr["NAXIS1"], buff_hdr["NAXIS2"]),
                                           dtype=np.float32)
                    for ____ in range(len(sel_table_4)):
                        dark_images[____] = read_fits_file(sel_table_4["FILEPATH"][____])[1]

                    mdark_data = np.nanmean(sigma_clip(dark_images, sigma=3, maxiters=5, masked=False, axis=0), axis=0)
                    buff_hdr["DATE-OBS"] = ""
                    buff_hdr["ALPHA"] = ""
                    buff_hdr["DELTA"] = ""
                    buff_hdr["CCD-TEMP"] = (mean_ccd_temp, "mean temperature of sum of expositions [C]")
                    buff_hdr["IMAGETYP"] = "Master Dark"
                    mdark_filename = f"MasterDark_B={_}_F={__}_E={___:.0f}_T={mean_ccd_temp:.0f}.fits"

                    match self.settings["MASTER_BITPIX"]:
                        case 16:
                            # uint16
                            mdark_data = np.round(mdark_data).astype(np.uint16)
                            buff_hdr['BITPIX'] = 16
                        case -64:
                            # float32
                            buff_hdr['BITPIX'] = -64
                        case _:
                            raise Exception(f"Unknown BITPIX for Master Dark: {self.settings['MASTER_BITPIX']}")

                    mdark_hdu = fits.PrimaryHDU(mdark_data, buff_hdr)
                    mdark_hdul = fits.HDUList([mdark_hdu])
                    mdark_hdul.writeto(self.settings["OUT_DIR"] + mdark_filename, overwrite=True)

                    imgs_type.append("Master Dark")
                    imgs_binning.append(_)
                    imgs_filter.append(__)
                    imgs_exptime.append(___)
                    imgs_ccdtemp.append(mean_ccd_temp)
                    imgs_filepath.append(self.settings["OUT_DIR"] + mdark_filename)

        if len(imgs_filepath) > 0:
            mdark_table = Table([imgs_type, imgs_binning, imgs_filter, imgs_exptime, imgs_ccdtemp, imgs_filepath],
                                names=("IMAGETYP", "BINNING", "FILTER", "EXPTIME", "CCD-TEMP", "FILEPATH"))
            self.images_table = vstack([self.images_table, mdark_table])
        print(f"Master darks created: {len(imgs_filepath)} total\n")

    def make_master_flat(self):
        check_out_directory(self.settings["OUT_DIR"])

        imgs_type = []
        imgs_binning = []
        imgs_filter = []
        imgs_exptime = []
        imgs_ccdtemp = []
        imgs_filepath = []

        sel_table_0 = self.images_table[np.where(((self.images_table["IMAGETYP"] == "Flat") |
                                                  (self.images_table["IMAGETYP"] == "Flat Field")))]
        if len(sel_table_0) == 0:
            print("Can't create master flat frames - No flat images found\n")
            return

        bin_factors = unique(sel_table_0, "BINNING")["BINNING"]
        for _ in bin_factors:
            sel_table_1 = sel_table_0[np.where(sel_table_0["BINNING"] == _)]
            filters = unique(sel_table_1, "FILTER")["FILTER"]
            for __ in filters:
                sel_table_2 = sel_table_1[np.where(sel_table_1["FILTER"] == __)]
                exptimes = unique(sel_table_2, "EXPTIME")["EXPTIME"]
                for ___ in exptimes:
                    sel_table_3 = sel_table_2[np.where(sel_table_2["EXPTIME"] == ___)]
                    mean_ccd_temp = np.round(np.mean(sel_table_3["CCD-TEMP"]), 3)
                    sel_table_4 = sel_table_3[np.where(abs(sel_table_3["CCD-TEMP"] - mean_ccd_temp) <=
                                                       self.settings["EPSILON_TEMP"])]
                    if len(sel_table_4) < self.settings["MIN_IMAGES_FOR_MASTER"]:
                        warn(f"Less than {self.settings['MIN_IMAGES_FOR_MASTER']} flat images with BINNING {_}, "
                             f"FILTER {__}, EXPTIME {___} and CCD-TEMP near {mean_ccd_temp} - skipped")
                        continue
                    print(f"Found {len(sel_table_4)} flat images with BINNING {_}, "
                          f"FILTER {__}, EXPTIME {___} and CCD-TEMP near {mean_ccd_temp}")
                    mdark_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Dark") &
                                                             (self.images_table["BINNING"] == _) &
                                                             (self.images_table["FILTER"] == __) &
                                                             (self.images_table["EXPTIME"] == ___) &
                                                             (abs(self.images_table["CCD-TEMP"] - mean_ccd_temp) <=
                                                              self.settings["EPSILON_TEMP"]))]
                    if len(mdark_table) == 0:
                        warn(f"No master dark images with BINNING {_}, FILTER {__}, EXPTIME {___} and"
                             f"CCD-TEMP near {mean_ccd_temp} - skipped")
                        continue

                    mdark_data = read_fits_file(mdark_table["FILEPATH"][0])[1]
                    buff_hdr = read_fits_file(sel_table_4["FILEPATH"][0])[0]
                    flat_images = np.zeros((len(sel_table_4), buff_hdr["NAXIS1"], buff_hdr["NAXIS2"]),
                                           dtype=np.float32)
                    for ____ in range(len(sel_table_4)):
                        flat_images[____] = read_fits_file(sel_table_4["FILEPATH"][___])[1]
                        flat_images[____] -= mdark_data
                        flat_images[____] /= np.mean(flat_images[___])

                    mflat_data = np.nanmean(sigma_clip(flat_images, sigma=3, maxiters=5, masked=False, axis=0), axis=0)
                    buff_hdr["DATE-OBS"] = ""
                    buff_hdr["ALPHA"] = ""
                    buff_hdr["DELTA"] = ""
                    buff_hdr["CCD-TEMP"] = (mean_ccd_temp, "mean temperature of sum of expositions")
                    buff_hdr["IMAGETYP"] = "Master Flat"
                    mflat_filename = f"MasterFlat_B={_}_F={__}_E={___:.0f}_T={mean_ccd_temp}.fits"

                    match self.settings["MASTER_BITPIX"]:
                        case 16:
                            # uint16
                            mflat_data = np.round(mflat_data).astype(np.uint16)
                            buff_hdr['BITPIX'] = 16
                        case -64:
                            # float32
                            buff_hdr['BITPIX'] = -64
                        case _:
                            raise Exception(f"Unknown BITPIX for Master Flat: {self.settings['MASTER_BITPIX']}")

                    buff_hdu = fits.PrimaryHDU(mflat_data, buff_hdr)
                    buff_hdul = fits.HDUList([buff_hdu])
                    buff_hdul.writeto(self.settings["OUT_DIR"] + mflat_filename, overwrite=True)

                    imgs_type.append("Master Flat")
                    imgs_binning.append(_)
                    imgs_filter.append(__)
                    imgs_exptime.append(___)
                    imgs_ccdtemp.append(mean_ccd_temp)
                    imgs_filepath.append(self.settings["OUT_DIR"] + mflat_filename)

        if len(imgs_filepath) > 0:
            mflat_table = Table([imgs_type, imgs_binning, imgs_filter, imgs_exptime, imgs_ccdtemp, imgs_filepath],
                                names=("IMAGETYP", "BINNING", "FILTER", "EXPTIME", "CCD-TEMP", "FILEPATH"))
            self.images_table = vstack([self.images_table, mflat_table])
        print(f"Master flats created: {len(imgs_filepath)} total\n")

    def apply_calibration(self):
        check_out_directory(self.settings["OUT_DIR"])

        __ = 0
        for _ in range(len(self.images_table)):
            if self.images_table["IMAGETYP"][_] in ["Object", "Light", "Light Frame"]:
                header, data = read_fits_file(self.images_table["FILEPATH"][_])

                match self.settings["CALIBRATION_TYPE"]:
                    case "DO_DARK":
                        # 'FILTER' check can be omitted
                        mdark_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Dark") &
                                                                 (self.images_table["BINNING"] ==
                                                                  self.images_table["BINNING"][_]) &
                                                                 (self.images_table["FILTER"] ==
                                                                  self.images_table["FILTER"][_]) &
                                                                 (self.images_table["EXPTIME"] ==
                                                                  self.images_table["EXPTIME"][_]) &
                                                                 (abs(self.images_table["CCD-TEMP"] -
                                                                      self.images_table["CCD-TEMP"][_]) <=
                                                                  self.settings["EPSILON_TEMP"]))]
                        if len(mdark_table) == 0:
                            warn(f"No fitting master dark images for {self.images_table['FILEPATH'][_]} - skipped")
                            continue
                        mdark_data = read_fits_file(mdark_table["FILEPATH"][0])[1]
                        print(f"File {self.images_table['FILEPATH'][_]}: Master Dark {mdark_table['FILEPATH'][0]}")
                        data = data.astype(np.float32) - mdark_data

                    case "DO_FLAT":
                        mflat_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Flat") &
                                                                 (self.images_table["BINNING"] ==
                                                                  self.images_table["BINNING"][_]) &
                                                                 (self.images_table["FILTER"] ==
                                                                  self.images_table["FILTER"][_]))]
                        if len(mflat_table) == 0:
                            warn(f"No fitting master flat images for {self.images_table['FILEPATH'][_]} - skipped")
                            continue
                        mflat_data = read_fits_file(mflat_table["FILEPATH"][0])[1]
                        print(f"File {self.images_table['FILEPATH'][_]}: Master Flat {mflat_table['FILEPATH'][0]}")
                        data = data.astype(np.float32) / mflat_data

                    case "DO_BOTH":
                        mdark_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Dark") &
                                                                 (self.images_table["BINNING"] ==
                                                                  self.images_table["BINNING"][_]) &
                                                                 (self.images_table["FILTER"] ==
                                                                  self.images_table["FILTER"][_]) &
                                                                 (self.images_table["EXPTIME"] ==
                                                                  self.images_table["EXPTIME"][_]) &
                                                                 (abs(self.images_table["CCD-TEMP"] -
                                                                      self.images_table["CCD-TEMP"][_]) <=
                                                                  self.settings["EPSILON_TEMP"]))]
                        if len(mdark_table) == 0:
                            warn(f"No fitting master dark images for {self.images_table['FILEPATH'][_]} - skipped")
                            continue
                        mflat_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Flat") &
                                                                 (self.images_table["BINNING"] ==
                                                                  self.images_table["BINNING"][_]) &
                                                                 (self.images_table["FILTER"] ==
                                                                  self.images_table["FILTER"][_]))]
                        if len(mflat_table) == 0:
                            warn(f"No fitting master flat images for {self.images_table['FILEPATH'][_]} - skipped")
                            continue
                        mdark_data = read_fits_file(mdark_table["FILEPATH"][0])[1]
                        mflat_data = read_fits_file(mflat_table["FILEPATH"][0])[1]
                        print(f"File {self.images_table['FILEPATH'][_]}: Master Dark {mdark_table['FILEPATH'][0]}, "
                              f"Master Flat {mflat_table['FILEPATH'][0]}")
                        data = data.astype(np.float32) - mdark_data
                        data /= mflat_data

                    case _:
                        raise Exception(f"Unknown CALIBRATION_TYPE value: {self.settings['CALIBRATION_TYPE']}")
                __ += 1

                match self.settings["CALIBRATED_BITPIX"]:
                    case 16:
                        # uint16
                        data = np.round(data).astype(np.uint16)
                        header['BITPIX'] = 16
                    case -64:
                        # float32
                        header['BITPIX'] = -64
                    case _:
                        raise Exception(f"Unknown BITPIX value for calibrated images: "
                                        f"{self.settings['CALIBRATED_BITPIX']}")
                header["HISTORY"] = "CALIBRATED"

                buff_hdu = fits.PrimaryHDU(data, header)
                buff_hdul = fits.HDUList([buff_hdu])

                # Slice off original file name
                buff_filepath = (self.settings["OUT_DIR"] +
                                 self.images_table["FILEPATH"][_][self.images_table["FILEPATH"][_].rfind("\\")+1:])
                buff_hdul.writeto(buff_filepath, overwrite=True)

        print(f"Images calibrated - {self.settings['CALIBRATION_TYPE']} mode, {__} total\n")


"""
Для дарков по идее не важен фильтр, но с его помощью можно учитывать
особенности генерации темнового тока у конкретной камеры

Для флэтов по идее не важна экспозиция, но этот лайфхак трудно учитывать
"""
