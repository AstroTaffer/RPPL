from datetime import datetime
from warnings import warn

from _rppl_ios import _RPPLSubIOS
from _rppl_aff import _RPPLSubFFA
from _rppl_clb import _RPPLSubCLB


class RPPLData(_RPPLSubIOS, _RPPLSubFFA, _RPPLSubCLB):
    # noinspection PyMissingConstructor
    def __init__(self, head_config=None, calc_config=None, **kwargs):
        self.head_pars = {"CONFIG_FILENAME": head_config,
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

        self.calc_pars = {"CONFIG_FILENAME": calc_config,
                          "IMAGES_DIRECTORY": None,
                          "IMAGES_BACKUP": True,
                          "FORCED_BACKUP": False,
                          "DELETE_NON-CALIBRATABLE": True,
                          "DELTA_T": None}

        self.images_table = None
        self.sdfimg_table = None

        if head_config is not None:
            self.read_head_config(f"{head_config}.json")
        if calc_config is not None:
            self.read_calc_config(f"{calc_config}.json")

        for _, __ in kwargs.items():
            if _ in self.head_pars:
                self.head_pars[_] = __
            elif _ in self.calc_pars:
                self.calc_pars[_] = __
            else:
                warn(f"Parameter {_} not recognised - skipped")

        self.calc_pars["IMAGES_DIRECTORY"] = self.calc_pars["IMAGES_DIRECTORY"].replace("\\", "/")
        if self.calc_pars["IMAGES_DIRECTORY"][-1] != "/":
            self.calc_pars["IMAGES_DIRECTORY"] += "/"

        for _ in self.head_pars.keys():
            if self.head_pars[_] is None:
                warn(f"Parameter {_} not set")
        for _ in self.calc_pars.keys():
            if self.calc_pars[_] is None:
                warn(f"Parameter {_} not set")

        print(f"{datetime.now()}\nNew RPPL data object created with following parameters:")
        print(self.head_pars)
        print(self.calc_pars)


"""

"""
