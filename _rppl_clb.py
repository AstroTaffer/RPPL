from os import remove
from warnings import warn

import numpy as np
import astropy.io.fits as fits
from astropy.table import unique, Table, vstack
from astropy.stats import sigma_clip

from utils import read_fits_file


class _RPPLSubCLB:
    def __init__(self):
        self.head_pars = None
        self.calc_pars = None
        self.images_table = None
        self.sdfimg_table = None

    def make_sdarks(self):
        imgs_filename = []
        imgs_filter = []
        imgs_type = []
        imgs_exptime = []
        imgs_ccdtemp = []
        imgs_sigma = []

        filters = unique(self.images_table, self.head_pars["FILTER"])[self.head_pars["FILTER"]]
        for _ in filters:
            sel_table = self.images_table[np.where((self.images_table[self.head_pars["FILTER"]] == _) &
                                                   ((self.images_table[self.head_pars["IMAGETYP"]] ==
                                                     self.head_pars["OBJECT"]) |
                                                    (self.images_table[self.head_pars["IMAGETYP"]] ==
                                                     self.head_pars["DARK"])))]
            exptimes = unique(sel_table, self.head_pars["EXPTIME"])[self.head_pars["EXPTIME"]]
            for __ in exptimes:
                sel_table = self.images_table[np.where((self.images_table[self.head_pars["IMAGETYP"]] ==
                                                        self.head_pars["DARK"]) &
                                                       (self.images_table[self.head_pars["FILTER"]] == _) &
                                                       (self.images_table[self.head_pars["EXPTIME"]] == __))]

                if len(sel_table) == 0:
                    warn(f"\nNo dark images with filter {_} and exptime {__} - skipped")
                    continue
                mean_ccd_temp = np.mean(sel_table[self.head_pars["CCD-TEMP"]])
                sel_table = sel_table[np.where(abs(sel_table[self.head_pars["CCD-TEMP"]] - mean_ccd_temp) <=
                                               self.calc_pars["DELTA_T"])]
                if len(sel_table) < 5:
                    warn(f"\nLess than 5 dark images with filter {_}, exptime {__} and"
                         f"temperature near {mean_ccd_temp} - skipped")
                    continue
                else:
                    print(f"\nTotal of {len(sel_table)} dark images with filter {_}, exptime {__} "
                          f"with mean temperature {mean_ccd_temp:.2f}")

                buff_hdr = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + sel_table["FILENAME"][0])[0]
                dark_images = np.zeros((len(sel_table), buff_hdr[self.head_pars["NAXIS1"]],
                                        buff_hdr[self.head_pars["NAXIS2"]]), dtype=np.float32)
                for ___ in range(len(sel_table)):
                    dark_images[___] = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] +
                                                      sel_table["FILENAME"][___])[1]

                sdark_image = np.nanmean(sigma_clip(dark_images, sigma=3, maxiters=5, masked=False, axis=0), axis=0)
                # FIXME: float32 or uint16 for sdark?
                buff_hdr[self.head_pars["DATE-OBS"]] = ""
                # FIXME: ADD MEAN DATE-OBS CALCULATION
                buff_hdr[self.head_pars["CCD-TEMP"]] = (mean_ccd_temp, "mean temperature of sum of expositions")
                buff_hdr[self.head_pars["IMAGETYP"]] = self.head_pars["SDARK"]

                sdark_filename = f"SDARK_{_}_{__}.fits"
                buff_hdu = fits.PrimaryHDU(sdark_image, buff_hdr)
                buff_hdul = fits.HDUList([buff_hdu])
                buff_hdul.writeto(self.calc_pars["IMAGES_DIRECTORY"] + sdark_filename, overwrite=True)

                imgs_filename.append(sdark_filename)
                imgs_filter.append(_)
                imgs_type.append(self.head_pars["SDARK"])
                imgs_exptime.append(__)
                imgs_ccdtemp.append(mean_ccd_temp)
                imgs_sigma.append(np.std(sdark_image))
        
        self.sdfimg_table = Table([imgs_filename, imgs_filter, imgs_type, imgs_exptime, imgs_ccdtemp, imgs_sigma],
                                  names=("FILENAME",
                                         self.head_pars["FILTER"],
                                         self.head_pars["IMAGETYP"],
                                         self.head_pars["EXPTIME"],
                                         self.head_pars["CCD-TEMP"],
                                         "SIGMA"))
        self.sdfimg_table.add_column(np.arange(1, len(self.sdfimg_table) + 1), index=0, name="ID")
        print(f"\nTotal of {len(self.sdfimg_table)} superdark images created:")
        print(self.sdfimg_table)

    def make_sflats(self):
        imgs_filename = []
        imgs_filter = []
        imgs_type = []
        imgs_exptime = []
        imgs_ccdtemp = []
        imgs_sigma = []

        filters = unique(self.images_table, self.head_pars["FILTER"])[self.head_pars["FILTER"]]
        for _ in filters:
            sel_table = self.images_table[np.where((self.images_table[self.head_pars["FILTER"]] == _) &
                                                   ((self.images_table[self.head_pars["IMAGETYP"]] ==
                                                     self.head_pars["OBJECT"]) |
                                                    (self.images_table[self.head_pars["IMAGETYP"]] ==
                                                     self.head_pars["FLAT"])))]
            exptimes = unique(sel_table, self.head_pars["EXPTIME"])[self.head_pars["EXPTIME"]]
            for __ in exptimes:
                sel_table = self.images_table[np.where((self.images_table[self.head_pars["IMAGETYP"]] ==
                                                        self.head_pars["FLAT"]) &
                                                       (self.images_table[self.head_pars["FILTER"]] == _) &
                                                       (self.images_table[self.head_pars["EXPTIME"]] == __))]

                if len(sel_table) == 0:
                    warn(f"\nNo flat images with filter {_} and exptime {__} - skipped")
                    continue
                mean_ccd_temp = np.mean(sel_table[self.head_pars["CCD-TEMP"]])
                sel_table = sel_table[np.where(abs(sel_table[self.head_pars["CCD-TEMP"]] - mean_ccd_temp) <=
                                               self.calc_pars["DELTA_T"])]
                if len(sel_table) < 5:
                    warn(f"\nLess than 5 flat images with filter {_}, exptime {__} and"
                         f"temperature near {mean_ccd_temp} - skipped")
                    continue
                else:
                    print(f"\nTotal of {len(sel_table)} flat images with filter {_}, exptime {__} "
                          f"with mean temperature {mean_ccd_temp:.2f}")

                sdark_table = self.sdfimg_table[np.where((self.sdfimg_table[self.head_pars["IMAGETYP"]] ==
                                                          self.head_pars["SDARK"]) &
                                                         (self.sdfimg_table[self.head_pars["EXPTIME"]] == __) &
                                                         (abs(self.sdfimg_table[self.head_pars["CCD-TEMP"]] -
                                                              mean_ccd_temp) <= self.calc_pars["DELTA_T"]))]
                if len(sdark_table) == 0:
                    warn(f"No suitable superdark image with filter {_}, exptime {__} and"
                         f"temperature near {mean_ccd_temp} - skipped")
                    continue

                # FIXME: ADD SELECTION BY NEAREST JD
                sdark_data = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + sdark_table["FILENAME"][0])[1]
                buff_hdr = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + sel_table["FILENAME"][0])[0]
                flat_images = np.zeros((len(sel_table), buff_hdr[self.head_pars["NAXIS1"]],
                                        buff_hdr[self.head_pars["NAXIS2"]]), dtype=np.float32)
                for ___ in range(len(sel_table)):
                    flat_images[___] = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] +
                                                      sel_table["FILENAME"][___])[1]
                    flat_images[___] -= sdark_data
                    flat_images[___] /= np.mean(flat_images[___])

                if self.calc_pars["AVERAGING_FUNC"] == "SIGMA-CLIP":
                    sflat_image = np.nanmean(sigma_clip(flat_images, sigma=3, maxiters=5, masked=False, axis=0), axis=0)
                elif self.calc_pars["AVERAGING_FUNC"] == "MEAN":
                    sflat_image = np.mean(flat_images, axis=0)
                elif self.calc_pars["AVERAGING_FUNC"] == "MEDIAN":
                    sflat_image = np.median(flat_images, axis=0)
                else:
                    warn("Unknown averaging function - skipped")
                    exit()

                buff_hdr[self.head_pars["DATE-OBS"]] = ""
                buff_hdr[self.head_pars["CCD-TEMP"]] = (mean_ccd_temp, "mean temperature of sum of expositions")
                buff_hdr[self.head_pars["IMAGETYP"]] = self.head_pars["SFLAT"]

                sflat_filename = f"SFLAT_{_}_{__}.fits"
                # noinspection PyUnboundLocalVariable
                # FIXME: Bad code - raise exception instead
                buff_hdu = fits.PrimaryHDU(sflat_image, buff_hdr)
                buff_hdul = fits.HDUList([buff_hdu])
                buff_hdul.writeto(self.calc_pars["IMAGES_DIRECTORY"] + sflat_filename, overwrite=True)

                imgs_filename.append(sflat_filename)
                imgs_filter.append(_)
                imgs_type.append(self.head_pars["SFLAT"])
                imgs_exptime.append(__)
                imgs_ccdtemp.append(mean_ccd_temp)
                imgs_sigma.append(np.std(sflat_image))

        sflat_table = Table([imgs_filename, imgs_filter, imgs_type, imgs_exptime, imgs_ccdtemp, imgs_sigma],
                            names=("FILENAME",
                                   self.head_pars["FILTER"],
                                   self.head_pars["IMAGETYP"],
                                   self.head_pars["EXPTIME"],
                                   self.head_pars["CCD-TEMP"],
                                   "SIGMA"))
        sflat_table.add_column(np.arange(len(self.sdfimg_table) + 1, len(self.sdfimg_table) + len(sflat_table) + 1),
                               index=0, name="ID")
        self.sdfimg_table = vstack([self.sdfimg_table, sflat_table])
        # FIXME: Somehow handle empty sflat_table
        print(f"\nTotal of {len(sflat_table)} superflat images created:")
        print(sflat_table)

    def apply_calibration(self):
        _ = 0
        while _ < len(self.images_table):
            if self.images_table[self.head_pars["IMAGETYP"]][_] != self.head_pars["OBJECT"]:
                _ += 1
                continue
            header, data = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + self.images_table["FILENAME"][_])
            sdark_table = self.sdfimg_table[np.where((self.sdfimg_table[self.head_pars["IMAGETYP"]] ==
                                                      self.head_pars["SDARK"]) &
                                                     (self.sdfimg_table[self.head_pars["EXPTIME"]] ==
                                                      header[self.head_pars["EXPTIME"]]) &
                                                     (abs(self.sdfimg_table[self.head_pars["CCD-TEMP"]] -
                                                          header[self.head_pars["CCD-TEMP"]]) <=
                                                      self.calc_pars["DELTA_T"]))]
            sflat_table = self.sdfimg_table[np.where((self.sdfimg_table[self.head_pars["IMAGETYP"]] ==
                                                      self.head_pars["SFLAT"]) &
                                                     (self.sdfimg_table[self.head_pars["FILTER"]] ==
                                                      header[self.head_pars["FILTER"]]) &
                                                     (abs(self.sdfimg_table[self.head_pars["CCD-TEMP"]] -
                                                          header[self.head_pars["CCD-TEMP"]]) <=
                                                      self.calc_pars["DELTA_T"]))]
            if len(sdark_table) == 0:
                warn(f"No suitable superdark image for image {self.images_table['FILENAME'][_]} - excluded")
                self.images_table.remove_row(_)
                if self.calc_pars["DELETE_NON-CALIBRATABLE"]:
                    remove(self.calc_pars["IMAGES_DIRECTORY"] + self.images_table["FILENAME"][_])
                    print("Image removed")
                continue
            elif len(sflat_table) == 0:
                warn(f"No suitable superflat image for image {self.images_table['FILENAME'][_]} - excluded")
                self.images_table.remove_row(_)
                if self.calc_pars["DELETE_NON-CALIBRATABLE"]:
                    remove(self.calc_pars["IMAGES_DIRECTORY"] + self.images_table["FILENAME"][_])
                    print("Image removed")
                continue
            else:
                sdark_data = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + sdark_table["FILENAME"][0])[1]
                sflat_data = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + sflat_table["FILENAME"][0])[1]
                print(f"FILE #{self.images_table['ID'][_]} - {self.images_table['FILENAME'][_]}, "
                      f"SDARK - {sdark_table['FILENAME'][0]}, SFLAT - {sflat_table['FILENAME'][0]}")

            data = data.astype(np.float32) - sdark_data
            data /= sflat_data
            data = np.around(data).astype(np.uint16)

            header['BITPIX'] = 16
            header["HISTORY"] = 'CALIBRATED'
            buff_hdu = fits.PrimaryHDU(data, header)
            buff_hdul = fits.HDUList([buff_hdu])
            buff_hdul.writeto(self.calc_pars["IMAGES_DIRECTORY"] + self.images_table["FILENAME"][_], overwrite=True)
            _ += 1
        print(f"Total of {len(self.images_table)} calibrated")

    def calibrate_images(self):
        self.make_sdarks()
        self.make_sflats()
        self.apply_calibration()


"""
Даркам не важен фильтр, флетам не важна экспозиция -> возможен выбор при калибровке
Возможно, этот выбор лучше сделать по наибольшему числу кадров, задействованных в создании супер-снимков.
Можно упростить перебор по сочетанию фильтр/экспозиция

Когда заработает возврат строк по астропаевскому мультииндексу, его стоит внедрить
"""