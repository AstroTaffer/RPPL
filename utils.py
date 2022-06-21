import json

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np


# FIXME: This function requires manual updating with every new parameter introduced. That's actually a bad solution
# LAST UPDATE 21-06 12:00
def create_empty_configs(head_filename="new_head_config", calc_filename="new_calc_config"):
    head_pars = {"CONFIG_FILENAME": head_filename,
                 "FILTER": None,
                 "IMAGETYP": None,
                 "EXPTIME": None,
                 "CCD-TEMP": None,
                 "NAXIS1": None,
                 "NAXIS2": None,
                 "DATE-OBS": None,

                 "OBJECT": None,
                 "BIAS": None,
                 "DARK": None,
                 "FLAT": None,
                 "SDARK": None,
                 "SFLAT": None}

    calc_pars = {"CONFIG_FILENAME": calc_filename,
                 "IMAGES_DIRECTORY": None,
                 "IMAGES_BACKUP": True,
                 "FORCED_BACKUP": False,
                 "DELETE_NON-CALIBRATABLE": True,
                 "DELTA_T": None,
                 "AVERAGING_FUNC": "SIGMA-CLIP"}

    with open(f"{head_filename}.json", "w") as head_confile:
        json.dump(head_pars, head_confile, indent=4)
    with open(f"{calc_filename}.json", "w") as calc_confile:
        json.dump(calc_pars, calc_confile, indent=4)


def read_fits_file(file_path):
    with fits.open(file_path) as hdu_list:
        hdu_list.verify("fix")
        img_header = hdu_list[0].header
        img_data = hdu_list[0].data
    return img_header, img_data


def plot_results():
    labels = ['Mean', 'Median', 'Sigma-clip']

    sec1 = [0.017447633668780327, 0.017519615590572357, 0.017447983492020033]
    sec2 = [0.016857948154211044,  0.01734243705868721, 0.016479508754926216]

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(dpi=150)
    rects1 = ax.bar(x - width / 2, sec1, width, label="EXPTIME = 1 sec")
    rects2 = ax.bar(x + width / 2, sec2, width, label="EXPTIME = 2 sec")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(r"$\sigma$")
    ax.set_title("r filter")
    ax.set_xticks(x, labels)
    ax.legend()

    ax.bar_label(rects1, padding=3, fmt="%.5f")
    ax.bar_label(rects2, padding=3, fmt="%.5f")

    fig.tight_layout()

    plt.ylim(0, 0.025)

    plt.savefig(f"r_res.png")
    plt.close(fig)


"""
r - sigma_clip
 ID     FILENAME     FILTER ... EXPTIME      CCD-TEMP             SIGMA        
--- ---------------- ------ ... ------- ------------------ --------------------
  3 SFLAT_r_1.0.fits      r ...     1.0           -50.0125 0.017447983492020033
  4 SFLAT_r_2.0.fits      r ...     2.0 -49.99431818181818 0.016479508754926216
       SIGMA        
--------------------


"""
