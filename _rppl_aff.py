from os import listdir, mkdir
from os.path import isdir
from shutil import copy2

import numpy as np
from astropy.table import Table

from utils import read_fits_file


class _RPPLSubFFA:
    def __init__(self):
        self.head_pars = None
        self.calc_pars = None
        self.images_table = None

    def account_fits_files(self):
        dir_content = listdir(self.calc_pars["IMAGES_DIRECTORY"])
        _ = 0
        while _ < len(dir_content):
            image_name = dir_content[_]
            if not (image_name.count(".fits") or image_name.count(".fit") or image_name.count(".fts")):
                dir_content.remove(image_name)
            else:
                _ += 1

        if self.calc_pars["IMAGES_BACKUP"]:
            self._backup_fits_files(dir_content)

        imgs_filename = []
        imgs_filter = []
        imgs_type = []
        imgs_exptime = []
        imgs_ccdtemp = []
        for _ in dir_content:
            header = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + _)[0]
            if not(header[self.head_pars["IMAGETYP"]] in [self.head_pars["BIAS"],
                                                          self.head_pars["DARK"],
                                                          self.head_pars["FLAT"],
                                                          self.head_pars["OBJECT"]]):
                # must be SDark or SFlat - none of them are allowed here
                continue
            imgs_filename.append(_)
            imgs_filter.append(header[self.head_pars["FILTER"]])
            # FIXME: Bias and Dark files sometimes have no 'FILTER' field
            # FIXME: I should test how astropy.io.fits handles missing header keys
            imgs_type.append(header[self.head_pars["IMAGETYP"]])
            imgs_exptime.append(header[self.head_pars["EXPTIME"]])
            imgs_ccdtemp.append(header[self.head_pars["CCD-TEMP"]])
        imgs_id = np.arange(1, len(imgs_filename) + 1)
        self.images_table = Table([imgs_id, imgs_filename, imgs_filter, imgs_type, imgs_exptime, imgs_ccdtemp],
                                  names=("ID", "FILENAME",
                                         self.head_pars["FILTER"],
                                         self.head_pars["IMAGETYP"],
                                         self.head_pars["EXPTIME"],
                                         self.head_pars["CCD-TEMP"]))

        print(f"\nTotal of {len(dir_content)} images accounted:")
        print(self.images_table)

    def _backup_fits_files(self, files_list):
        bkp_path = self.calc_pars["IMAGES_DIRECTORY"] + "backup/"
        if not isdir(bkp_path):
            mkdir(bkp_path)

        # Warning! Some metadata may be lost due to Python inbuilt file copying functions limitations.
        bkp_counter = 0
        bkp_content = listdir(bkp_path)
        for _ in files_list:
            if _ not in bkp_content or self.calc_pars["FORCED_BACKUP"]:
                # FIXME: copy2 can't parse cyrillic symbols
                copy2(self.calc_pars["IMAGES_DIRECTORY"] + _, bkp_path + _)
                bkp_counter += 1
        print(f"\nBackup procedured - {bkp_counter} new files created, {len(files_list) - bkp_counter} already existed")


"""

"""
