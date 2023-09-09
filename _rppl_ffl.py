from os import listdir

from astropy.table import Table

from utils import read_fits_file


class _Lister:
    def __init__(self):
        self.settings = None
        self.images_table = None

    def list_fits_files(self):
        imgs_type = []
        imgs_binning = []
        imgs_filter = []
        imgs_exptime = []
        imgs_ccdtemp = []
        imgs_filepath = []

        for _ in self.settings["IN_DIRS"]:
            for __ in listdir(_):
                if __.count(".fits") or __.count(".fit") or __.count(".fts"):
                    header = read_fits_file(_ + __)[0]
                    imgs_type.append(header["IMAGETYP"])
                    imgs_binning.append(header["BINNING"])
                    try:
                        imgs_filter.append(header["FILTER"])
                    except KeyError:
                        imgs_filter.append("UNKNW")
                    imgs_exptime.append(header["EXPTIME"])
                    imgs_ccdtemp.append(round(header["CCD-TEMP"], 3))
                    imgs_filepath.append(_ + __)

        self.images_table = Table([imgs_type, imgs_binning, imgs_filter, imgs_exptime, imgs_ccdtemp, imgs_filepath],
                                  names=("IMAGETYP", "BINNING", "FILTER", "EXPTIME", "CCD-TEMP", "FILEPATH"))
        print(f"FITS files listed: {len(self.images_table)} found\n")
