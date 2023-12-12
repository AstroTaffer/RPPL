import json
from os.path import exists
import os
import subprocess
import numpy as np
from astropy.io import ascii
from sqlalchemy import create_engine
from sqlalchemy import URL
import astropy.io.fits as fits


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


def do_sex(input_file):
    # cwd = 'F:\\'
    cwd = os.getcwd() + '\\'
    Sex = cwd + 'Sex\Extract.exe '
    dSex = ' -c ' + cwd + 'Sex\default.sex'
    dPar = ' -PARAMETERS_NAME ' + cwd + 'Sex\default.par'
    dFilt = ' -FILTER_NAME ' + cwd + r'Sex\tophat_2.5_3x3.conv'
    NNW = ' -STARNNW_NAME ' + cwd + 'Sex\default.nnw'
    output_file = input_file.replace('fits.gz', 'cat')

    shell = Sex + "\"" + input_file + "\"" + dSex + dPar + dFilt + NNW + ' -CATALOG_NAME ' + "\"" + output_file + "\""
    print(shell)
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    child = subprocess.run(shell, timeout=60, startupinfo=startupinfo)

    if child.returncode == 0 and os.path.isfile(output_file):
        print('Ok')
    else:
        print('Error')
        return 0, 0, 0, 0, ''
    tbl = ascii.read(output_file)
    med_fwhm = np.round(np.median(tbl['FWHM_IMAGE']), 2)
    med_ell = np.round(np.median(tbl['ELLIPTICITY']), 2)
    med_bkg = np.round(np.median(tbl['BACKGROUND']), 2)
    med_zeropoi = np.round(np.median(tbl['ZEROPOI']), 2)
    return med_fwhm, med_ell, med_bkg, med_zeropoi, output_file


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