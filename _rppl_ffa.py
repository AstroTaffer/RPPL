from os import listdir, mkdir
from os.path import isdir
from shutil import copy2
from datetime import datetime

import numpy as np
from astropy.table import Table

from utils import read_fits_file


class _RPPLSubFFA:
    def __init__(self):
        self.head_pars = {}
        self.calc_pars = {}
        self.images_table = None

    def fits_files_accounting(self):
        dir_content = listdir(self.calc_pars["IMAGES_DIRECTORY"])
        _ = 0
        while _ < len(dir_content):
            image_name = dir_content[_]
            if not (image_name.count(".fits") or image_name.count(".fit") or image_name.count(".fts")):
                dir_content.remove(image_name)
            _ += 1

        if self.calc_pars["IMAGES_BACKUP"]:
            self._fits_files_backup(dir_content)

        imgs_id = np.arange(1, len(dir_content) + 1)
        imgs_filename = []
        imgs_filter = []
        imgs_type = []
        imgs_exptime = []
        imgs_ccdtemp = []
        for _ in dir_content:
            header = read_fits_file(self.calc_pars["IMAGES_DIRECTORY"] + _)[0]
            imgs_filename.append(_)
            imgs_filter.append(header[self.head_pars["FILTER"]])
            # FIXME Bias and Dark files sometimes have no 'FILTER' field
            # FIXME: I should test how astropy.io.fits handles missing header keys
            imgs_type.append(header[self.head_pars["IMAGETYP"]])
            imgs_exptime.append(header[self.head_pars["EXPTIME"]])
            imgs_ccdtemp.append(header[self.head_pars["CCD-TEMP"]])
        self.images_table = Table([imgs_id, imgs_filename, imgs_filter, imgs_type, imgs_exptime, imgs_ccdtemp],
                                  names=("ID", "FILENAME",
                                         self.head_pars["FILTER"],
                                         self.head_pars["IMAGETYP"],
                                         self.head_pars["EXPTIME"],
                                         self.head_pars["CCD-TEMP"]))

        print(f"{datetime.now()}\nTotal of {len(dir_content)} images accounted:")
        print(self.images_table)

    def _fits_files_backup(self, files_list):
        bkp_path = self.calc_pars["IMAGES_DIRECTORY"] + "backup/"
        if not isdir(bkp_path):
            mkdir("backup/")

        # FIXME: REMOVE TEMPO_COUNTER
        # FIXME: IF COPY2 IS TOO SLOW, CHOOSE ANOTHER METHOD
        # FIXME: WARN USER ABOUT LOSING METADATA
        tmp_counter = 1
        bkp_content = listdir(bkp_path)
        for _ in files_list:
            if _ not in bkp_content or self.calc_pars["FORCED_BACKUP"]:
                copy2(self.calc_pars["IMAGES_DIRECTORY"] + _, bkp_path + _)
                print(f"FIle #{tmp_counter} {_} copied")
            else:
                print(f"File #{tmp_counter} {_} already exists")
            tmp_counter += 1


"""

"""
