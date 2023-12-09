import numpy as np
import astropy.io.fits as fits
from astropy.stats import sigma_clip

from utils import read_fits_file


def apply_dark(raw_frame_fp, mdark_frame_fp, out_frame_fp):
    try:
        raw_frame_header, raw_frame_data = read_fits_file(raw_frame_fp)
        mdark_frame_data = read_fits_file(mdark_frame_fp)[1]

        raw_frame_data = raw_frame_data.astype(np.float32) - mdark_frame_data

        raw_frame_data = np.round(raw_frame_data).astype(np.uint16)
        raw_frame_header["BITPIX"] = (16, "bits per data value")
        raw_frame_header["HISTORY"] = "CALIBRATED"

        out_frame_hdu = fits.PrimaryHDU(raw_frame_data, raw_frame_header)
        out_frame_hdul = fits.HDUList([out_frame_hdu])
        out_frame_hdul.writeto(out_frame_fp, overwrite=True)
        return True
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return False


def apply_flat(raw_frame_fp, mflat_frame_fp, out_frame_fp):
    try:
        raw_frame_header, raw_frame_data = read_fits_file(raw_frame_fp)
        mflat_frame_data = read_fits_file(mflat_frame_fp)[1]

        raw_frame_data = raw_frame_data.astype(np.float32) / mflat_frame_data

        raw_frame_data = np.round(raw_frame_data).astype(np.uint16)
        raw_frame_header["BITPIX"] = (16, "bits per data value")
        raw_frame_header["HISTORY"] = "CALIBRATED"

        out_frame_hdu = fits.PrimaryHDU(raw_frame_data, raw_frame_header)
        out_frame_hdul = fits.HDUList([out_frame_hdu])
        out_frame_hdul.writeto(out_frame_fp, overwrite=True)
        return True
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return False


def make_master_dark(dark_frames_fp, mean_ccd_temp, creation_date, out_frame_fp):
    try:
        buff_frame_hdr = read_fits_file(dark_frames_fp[0])[0]
        dark_frames_data = np.zeros((len(dark_frames_fp), buff_frame_hdr["NAXIS1"], buff_frame_hdr["NAXIS2"]),
                                    dtype=np.float32)

        for _ in range(len(dark_frames_fp)):
            dark_frames_data[_] = read_fits_file(dark_frames_fp[_])[1]

        mdark_frame_data = np.nanmean(sigma_clip(dark_frames_data, sigma=3, maxiters=5, masked=False, axis=0), axis=0)

        del buff_frame_hdr["DATE-OBS"]
        del buff_frame_hdr["ALPHA"]
        del buff_frame_hdr["DELTA"]
        buff_frame_hdr["BITPIX"] = (-64, "bits per data value")
        buff_frame_hdr["DATE"] = (creation_date, "date-time of file creation")
        buff_frame_hdr["CCD-TEMP"] = (mean_ccd_temp, "mean temperature of sum of expositions [C]")
        buff_frame_hdr["IMAGETYP"] = "Master Dark"

        mdark_frame_hdu = fits.PrimaryHDU(mdark_frame_data, buff_frame_hdr)
        mdark_frame_hdul = fits.HDUList([mdark_frame_hdu])
        mdark_frame_hdul.writeto(out_frame_fp, overwrite=True)
        return True
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return False


def make_master_flat(flat_frames_fp, mdark_frame_fp, mean_ccd_temp, creation_date, out_frame_fp):
    try:
        mdark_frame_data = read_fits_file(mdark_frame_fp)[1]

        buff_frame_hdr = read_fits_file(flat_frames_fp[0])[0]
        flat_frames_data = np.zeros((len(flat_frames_fp), buff_frame_hdr["NAXIS1"], buff_frame_hdr["NAXIS2"]),
                                    dtype=np.float32)

        for _ in range(len(flat_frames_fp)):
            flat_frames_data[_] = read_fits_file(flat_frames_fp[_])[1]
            flat_frames_data[_] -= mdark_frame_data
            flat_frames_data[_] /= np.mean(flat_frames_data[_])

        mflat_frame_data = np.nanmean(sigma_clip(flat_frames_data, sigma=3, maxiters=5, masked=False, axis=0), axis=0)

        del buff_frame_hdr["DATE-OBS"]
        del buff_frame_hdr["ALPHA"]
        del buff_frame_hdr["DELTA"]
        buff_frame_hdr["BITPIX"] = (-64, "bits per data value")
        buff_frame_hdr["DATE"] = (creation_date, "date-time of file creation")
        buff_frame_hdr["CCD-TEMP"] = (mean_ccd_temp, "mean temperature of sum of expositions [C]")
        buff_frame_hdr["IMAGETYP"] = "Master Flat"

        mflat_frame_hdu = fits.PrimaryHDU(mflat_frame_data, buff_frame_hdr)
        mflat_frame_hdul = fits.HDUList([mflat_frame_hdu])
        mflat_frame_hdul.writeto(out_frame_fp, overwrite=True)
        return True
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return False


# ----- KEEP A COPY OF OLD CODE JUST IN CASE ---- #

# class _Calibrator_Old:
#     def __init__(self):
#         self.settings = None
#         self.images_table = None
# 
#     def old_make_master_dark(self):
#         check_out_directory(self.settings["OUT_DIR"])
# 
#         imgs_type = []
#         imgs_xbinning = []
#         imgs_ybinning = []
#         imgs_filter = []
#         imgs_exptime = []
#         imgs_ccdtemp = []
#         imgs_filepath = []
# 
#         sel_table_0 = self.images_table[np.where(((self.images_table["IMAGETYP"] == "Dark") |
#                                                   (self.images_table["IMAGETYP"] == "Dark Frame")))]
#         if len(sel_table_0) == 0:
#             print("Can't create master dark frames - No dark images found\n")
#             return
# 
#         xbin_factors = unique(sel_table_0, "XBINNING")["XBINNING"]
#         for _ in xbin_factors:
#             sel_table_01 = sel_table_0[np.where(sel_table_0["XBINNING"] == _)]
#             ybin_factors = unique(sel_table_01, "YBINNING")["YBINNING"]
#             for __ in ybin_factors:
#                 sel_table_1 = sel_table_01[np.where(sel_table_01["YBINNING"] == __)]
#                 filters = unique(sel_table_1, "FILTER")["FILTER"]
#                 for ___ in filters:
#                     sel_table_2 = sel_table_1[np.where(sel_table_1["FILTER"] == ___)]
#                     exptimes = unique(sel_table_2, "EXPTIME")["EXPTIME"]
#                     for ____ in exptimes:
#                         sel_table_3 = sel_table_2[np.where(sel_table_2["EXPTIME"] == ____)]
#                         mean_ccd_temp = np.round(np.mean(sel_table_3["CCD-TEMP"]), 3)
#                         sel_table_4 = sel_table_3[np.where(abs(sel_table_3["CCD-TEMP"] - mean_ccd_temp) <=
#                                                            self.settings["EPSILON_TEMP"])]
#                         if len(sel_table_4) < self.settings["MIN_IMAGES_FOR_MASTER"]:
#                             warn(f"Less than {self.settings['MIN_IMAGES_FOR_MASTER']} dark images with XBINNING {_}, "
#                                  f"YBINNING {__}, FILTER {___}, EXPTIME {____} "
#                                  f"and CCD-TEMP near {mean_ccd_temp} - skipped")
#                             continue
#                         print(f"Found {len(sel_table_4)} dark images with XBINNING {_}, YBINNING {__}, "
#                               f"FILTER {___}, EXPTIME {____} and CCD-TEMP near {mean_ccd_temp}")
#     
#                         buff_hdr = read_fits_file(sel_table_4["FILEPATH"][0])[0]
#                         dark_images = np.zeros((len(sel_table_4), buff_hdr["NAXIS1"], buff_hdr["NAXIS2"]),
#                                                dtype=np.float32)
#                         for _____ in range(len(sel_table_4)):
#                             dark_images[_____] = read_fits_file(sel_table_4["FILEPATH"][_____])[1]
#     
#                         mdark_data = np.nanmean(sigma_clip(dark_images, sigma=3, maxiters=5, masked=False, axis=0),
#                                                 axis=0)
#                         buff_hdr['BITPIX'] = -64
#                         buff_hdr["DATE-OBS"] = ""
#                         buff_hdr["ALPHA"] = ""
#                         buff_hdr["DELTA"] = ""
#                         buff_hdr["CCD-TEMP"] = (mean_ccd_temp, "mean temperature of sum of expositions [C]")
#                         buff_hdr["IMAGETYP"] = "Master Dark"
#                         mdark_filename = f"MasterDark_XB={_}_YB={__}_F={___}_E={____:.0f}_T={mean_ccd_temp:.0f}.fits"
#     
#                         mdark_hdu = fits.PrimaryHDU(mdark_data, buff_hdr)
#                         mdark_hdul = fits.HDUList([mdark_hdu])
#                         mdark_hdul.writeto(self.settings["OUT_DIR"] + mdark_filename, overwrite=True)
#     
#                         imgs_type.append("Master Dark")
#                         imgs_xbinning.append(_)
#                         imgs_ybinning.append(__)
#                         imgs_filter.append(___)
#                         imgs_exptime.append(____)
#                         imgs_ccdtemp.append(mean_ccd_temp)
#                         imgs_filepath.append(self.settings["OUT_DIR"] + mdark_filename)
# 
#         if len(imgs_filepath) > 0:
#             mdark_table = Table([imgs_type, imgs_xbinning, imgs_ybinning, imgs_filter, imgs_exptime, imgs_ccdtemp,
#                                  imgs_filepath],
#                                 names=("IMAGETYP", "XBINNING", "YBINNING", "FILTER", "EXPTIME", "CCD-TEMP", "FILEPATH"))
#             self.images_table = vstack([self.images_table, mdark_table])
#         print(f"Master darks created: {len(imgs_filepath)} total\n")
# 
#     def old_make_master_flat(self):
#         check_out_directory(self.settings["OUT_DIR"])
# 
#         imgs_type = []
#         imgs_xbinning = []
#         imgs_ybinning = []
#         imgs_filter = []
#         imgs_exptime = []
#         imgs_ccdtemp = []
#         imgs_filepath = []
# 
#         sel_table_0 = self.images_table[np.where(((self.images_table["IMAGETYP"] == "Flat") |
#                                                   (self.images_table["IMAGETYP"] == "Flat Field")))]
#         if len(sel_table_0) == 0:
#             print("Can't create master flat frames - No flat images found\n")
#             return
# 
#         xbin_factors = unique(sel_table_0, "XBINNING")["XBINNING"]
#         for _ in xbin_factors:
#             sel_table_01 = sel_table_0[np.where(sel_table_0["XBINNING"] == _)]
#             ybin_factors = unique(sel_table_01, "YBINNING")["YBINNING"]
#             for __ in ybin_factors:
#                 sel_table_1 = sel_table_01[np.where(sel_table_01["YBINNING"] == __)]
#                 filters = unique(sel_table_1, "FILTER")["FILTER"]
#                 for ___ in filters:
#                     sel_table_2 = sel_table_1[np.where(sel_table_1["FILTER"] == ___)]
#                     exptimes = unique(sel_table_2, "EXPTIME")["EXPTIME"]
#                     for ____ in exptimes:
#                         sel_table_3 = sel_table_2[np.where(sel_table_2["EXPTIME"] == ____)]
#                         mean_ccd_temp = np.round(np.mean(sel_table_3["CCD-TEMP"]), 3)
#                         sel_table_4 = sel_table_3[np.where(abs(sel_table_3["CCD-TEMP"] - mean_ccd_temp) <=
#                                                            self.settings["EPSILON_TEMP"])]
#                         if len(sel_table_4) < self.settings["MIN_IMAGES_FOR_MASTER"]:
#                             warn(f"Less than {self.settings['MIN_IMAGES_FOR_MASTER']} flat images with XBINNING {_}, "
#                                  f"YBINNING {__}, FILTER {___}, EXPTIME {____} "
#                                  f"and CCD-TEMP near {mean_ccd_temp} - skipped")
#                             continue
#                         print(f"Found {len(sel_table_4)} flat images with XBINNING {_}, YBINNING {__}, "
#                               f"FILTER {___}, EXPTIME {____} and CCD-TEMP near {mean_ccd_temp}")
#                         mdark_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Dark") &
#                                                                  (self.images_table["XBINNING"] == _) &
#                                                                  (self.images_table["YBINNING"] == __) &
#                                                                  (self.images_table["FILTER"] == ___) &
#                                                                  (self.images_table["EXPTIME"] == ____) &
#                                                                  (abs(self.images_table["CCD-TEMP"] - mean_ccd_temp) <=
#                                                                   self.settings["EPSILON_TEMP"]))]
#                         if len(mdark_table) == 0:
#                             warn(f"No master dark images with XBINNING {_}, YBINNING {__}, FILTER {___}, "
#                                  f"EXPTIME {____} and CCD-TEMP near {mean_ccd_temp} - skipped")
#                             continue
# 
#                         mdark_data = read_fits_file(mdark_table["FILEPATH"][0])[1]
#                         buff_hdr = read_fits_file(sel_table_4["FILEPATH"][0])[0]
#                         flat_images = np.zeros((len(sel_table_4), buff_hdr["NAXIS1"], buff_hdr["NAXIS2"]),
#                                                dtype=np.float32)
#                         for _____ in range(len(sel_table_4)):
#                             flat_images[_____] = read_fits_file(sel_table_4["FILEPATH"][_____])[1]
#                             flat_images[_____] -= mdark_data
#                             flat_images[_____] /= np.mean(flat_images[_____])
# 
#                         mflat_data = np.nanmean(sigma_clip(flat_images, sigma=3, maxiters=5, masked=False, axis=0),
#                                                 axis=0)
#                         buff_hdr['BITPIX'] = -64
#                         buff_hdr["DATE-OBS"] = ""
#                         buff_hdr["ALPHA"] = ""
#                         buff_hdr["DELTA"] = ""
#                         buff_hdr["CCD-TEMP"] = (mean_ccd_temp, "mean temperature of sum of expositions")
#                         buff_hdr["IMAGETYP"] = "Master Flat"
#                         mflat_filename = f"MasterFlat_XB={_}_YB={__}_F={___}.fits"
# 
#                         buff_hdu = fits.PrimaryHDU(mflat_data, buff_hdr)
#                         buff_hdul = fits.HDUList([buff_hdu])
#                         buff_hdul.writeto(self.settings["OUT_DIR"] + mflat_filename, overwrite=True)
# 
#                         imgs_type.append("Master Flat")
#                         imgs_xbinning.append(_)
#                         imgs_ybinning.append(__)
#                         imgs_filter.append(___)
#                         imgs_exptime.append(____)
#                         imgs_ccdtemp.append(mean_ccd_temp)
#                         imgs_filepath.append(self.settings["OUT_DIR"] + mflat_filename)
# 
#         if len(imgs_filepath) > 0:
#             mflat_table = Table([imgs_type, imgs_xbinning, imgs_ybinning, imgs_filter, imgs_exptime, imgs_ccdtemp,
#                                  imgs_filepath],
#                                 names=("IMAGETYP", "XBINNING", "YBINNING", "FILTER", "EXPTIME", "CCD-TEMP", "FILEPATH"))
#             self.images_table = vstack([self.images_table, mflat_table])
#         print(f"Master flats created: {len(imgs_filepath)} total\n")
# 
#     def old_apply_calibration(self):
#         check_out_directory(self.settings["OUT_DIR"])
# 
#         __ = 0
#         for _ in range(len(self.images_table)):
#             if self.images_table["IMAGETYP"][_] in ["Object", "Light", "Light Frame"]:
#                 header, data = read_fits_file(self.images_table["FILEPATH"][_])
# 
#                 match self.settings["CALIBRATION_TYPE"]:
#                     case "DO_DARK":
#                         # 'FILTER' check can be omitted
#                         mdark_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Dark") &
#                                                                  (self.images_table["XBINNING"] ==
#                                                                   self.images_table["XBINNING"][_]) &
#                                                                  (self.images_table["YBINNING"] ==
#                                                                   self.images_table["YBINNING"][_]) &
#                                                                  (self.images_table["FILTER"] ==
#                                                                   self.images_table["FILTER"][_]) &
#                                                                  (self.images_table["EXPTIME"] ==
#                                                                   self.images_table["EXPTIME"][_]) &
#                                                                  (abs(self.images_table["CCD-TEMP"] -
#                                                                       self.images_table["CCD-TEMP"][_]) <=
#                                                                   self.settings["EPSILON_TEMP"]))]
#                         if len(mdark_table) == 0:
#                             warn(f"No fitting master dark images for {self.images_table['FILEPATH'][_]} - skipped")
#                             continue
#                         mdark_data = read_fits_file(mdark_table["FILEPATH"][0])[1]
#                         print(f"File {self.images_table['FILEPATH'][_]}: Master Dark {mdark_table['FILEPATH'][0]}")
#                         data = data.astype(np.float32) - mdark_data
# 
#                     case "DO_FLAT":
#                         mflat_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Flat") &
#                                                                  (self.images_table["XBINNING"] ==
#                                                                   self.images_table["XBINNING"][_]) &
#                                                                  (self.images_table["YBINNING"] ==
#                                                                   self.images_table["YBINNING"][_]) &
#                                                                  (self.images_table["FILTER"] ==
#                                                                   self.images_table["FILTER"][_]))]
#                         if len(mflat_table) == 0:
#                             warn(f"No fitting master flat images for {self.images_table['FILEPATH'][_]} - skipped")
#                             continue
#                         mflat_data = read_fits_file(mflat_table["FILEPATH"][0])[1]
#                         print(f"File {self.images_table['FILEPATH'][_]}: Master Flat {mflat_table['FILEPATH'][0]}")
#                         data = data.astype(np.float32) / mflat_data
# 
#                     case "DO_BOTH":
#                         mdark_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Dark") &
#                                                                  (self.images_table["XBINNING"] ==
#                                                                   self.images_table["XBINNING"][_]) &
#                                                                  (self.images_table["YBINNING"] ==
#                                                                   self.images_table["YBINNING"][_]) &
#                                                                  (self.images_table["FILTER"] ==
#                                                                   self.images_table["FILTER"][_]) &
#                                                                  (self.images_table["EXPTIME"] ==
#                                                                   self.images_table["EXPTIME"][_]) &
#                                                                  (abs(self.images_table["CCD-TEMP"] -
#                                                                       self.images_table["CCD-TEMP"][_]) <=
#                                                                   self.settings["EPSILON_TEMP"]))]
#                         if len(mdark_table) == 0:
#                             warn(f"No fitting master dark images for {self.images_table['FILEPATH'][_]} - skipped")
#                             continue
#                         mflat_table = self.images_table[np.where((self.images_table["IMAGETYP"] == "Master Flat") &
#                                                                  (self.images_table["XBINNING"] ==
#                                                                   self.images_table["XBINNING"][_]) &
#                                                                  (self.images_table["YBINNING"] ==
#                                                                   self.images_table["YBINNING"][_]) &
#                                                                  (self.images_table["FILTER"] ==
#                                                                   self.images_table["FILTER"][_]))]
#                         if len(mflat_table) == 0:
#                             warn(f"No fitting master flat images for {self.images_table['FILEPATH'][_]} - skipped")
#                             continue
#                         mdark_data = read_fits_file(mdark_table["FILEPATH"][0])[1]
#                         mflat_data = read_fits_file(mflat_table["FILEPATH"][0])[1]
#                         print(f"File {self.images_table['FILEPATH'][_]}: Master Dark {mdark_table['FILEPATH'][0]}, "
#                               f"Master Flat {mflat_table['FILEPATH'][0]}")
#                         data = data.astype(np.float32) - mdark_data
#                         data /= mflat_data
# 
#                     case _:
#                         raise Exception(f"Unknown CALIBRATION_TYPE value: {self.settings['CALIBRATION_TYPE']}")
#                 __ += 1
# 
#                 match self.settings["CALIBRATED_BITPIX"]:
#                     case 16:
#                         # uint16
#                         data = np.round(data).astype(np.uint16)
#                         header['BITPIX'] = 16
#                     case -64:
#                         # float32
#                         header['BITPIX'] = -64
#                     case _:
#                         raise Exception(f"Unknown BITPIX value for calibrated images: "
#                                         f"{self.settings['CALIBRATED_BITPIX']}")
#                 header["HISTORY"] = "CALIBRATED"
# 
#                 buff_hdu = fits.PrimaryHDU(data, header)
#                 buff_hdul = fits.HDUList([buff_hdu])
# 
#                 # Slice off original file name
#                 buff_filepath = (self.settings["OUT_DIR"] +
#                                  self.images_table["FILEPATH"][_][self.images_table["FILEPATH"][_].rfind("\\")+1:])
#                 buff_hdul.writeto(buff_filepath, overwrite=True)
# 
#         print(f"Images calibrated - {self.settings['CALIBRATION_TYPE']} mode, {__} total\n")


"""
Для дарков по идее не важен фильтр, но с его помощью можно учитывать
особенности генерации темнового тока у конкретной камеры

Для флэтов по идее не важна экспозиция, но этот лайфхак трудно учитывать
"""
