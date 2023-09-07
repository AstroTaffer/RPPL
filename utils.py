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


# Was made for quick and stupid fix --- please, do not use
def copy_header(ref_file, bad_file, new_dir):
    check_out_directory(new_dir)
    if not (ref_file.count(".fits") or ref_file.count(".fit") or ref_file.count(".fts")):
        return
    if not (bad_file.count(".fits") or bad_file.count(".fit") or bad_file.count(".fts")):
        return
    ref_head = read_fits_file(ref_file)[0]
    bad_head, bad_data = read_fits_file(bad_file)
    if ref_head["DATE-OBS"] != bad_head["DATE-OBS"]:
        return
    new_hdu = fits.PrimaryHDU(bad_data, ref_head)
    new_hdul = fits.HDUList([new_hdu])
    new_hdul.writeto(new_dir + bad_file[bad_file.rfind("\\")+1:], overwrite=True)
    print(f"{bad_file} fixed")
