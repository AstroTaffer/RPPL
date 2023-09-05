import json
from os.path import exists
from os import makedirs

import astropy.io.fits as fits


def read_fits_file(file_name):
    with fits.open(file_name) as hdu_list:
        hdu_list.verify("fix")
        header = hdu_list[0].header
        data = hdu_list[0].data
    return header, data


def check_out_directory(out_dir):
    if not exists(out_dir):
        makedirs(out_dir)


# LAST UPDATE 05-09-2023
def restore_default_config():
    settings = {
        "IN_DIRS": "",
        "CALIBRATION_TYPE": "DO_BOTH",
        "OUT_DIR": "",
        "EPSILON_TEMP": 1,
        "MIN_IMAGES_FOR_MASTER": 5,
        "MASTER_BITPIX": -64,
        "CALIBRATED_BITPIX": 16}
    with open(f"default_config.json", "w") as confile:
        json.dump(settings, confile, indent=4)
