import json

import astropy.io.fits as fits


# FIXME: This function requires manual updating with every new parameter introduced. That's actually a bad solution
# LAST UPDATE 30-05 09:50
def create_empty_configs(head_filename="new_head_config.json", calc_filename="new_calc_config.json"):
    head_pars = {"CONFIG_FILENAME": head_filename,
                 "FILTER": None,
                 "IMAGETYP": None,
                 "EXPTIME": None,
                 "CCD-TEMP": None,

                 "OBJECT": None,
                 "BIAS": None,
                 "DARK": None,
                 "FLAT": None}

    calc_pars = {"CONFIG_FILENAME": calc_filename,
                 "IMAGES_DIRECTORY": None,
                 "IMAGES_BACKUP": True,
                 "FORCED_BACKUP": False}

    with open(head_filename, "w") as head_confile:
        json.dump(head_pars, head_confile, indent=4)
    with open(calc_filename, "w") as calc_confile:
        json.dump(calc_pars, calc_confile, indent=4)


def read_fits_file(file_path):
    with fits.open(file_path) as hdu_list:
        hdu_list.verify("fix")
        img_header = hdu_list[0].header
        img_data = hdu_list[0].data
    return img_header, img_data


"""

"""
