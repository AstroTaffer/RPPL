from datetime import datetime
from warnings import warn

from _rppl_ios import _RPPLSubIOS
from _rppl_ffa import _RPPLSubFFA


class RPPLData(_RPPLSubIOS, _RPPLSubFFA):
    # noinspection PyMissingConstructor
    def __init__(self, head_config=None, calc_config=None, **kwargs):
        self.head_pars = {"CONFIG_FILENAME": head_config,
                          "FILTER": None,
                          "IMAGETYP": None,
                          "EXPTIME": None,
                          "CCD-TEMP": None,

                          "OBJECT": None,
                          "BIAS": None,
                          "DARK": None,
                          "FLAT": None}

        self.calc_pars = {"CONFIG_FILENAME": calc_config,
                          "IMAGES_DIRECTORY": None,
                          "IMAGES_BACKUP": True,
                          "FORCED_BACKUP": False}

        self.images_table = None

        if head_config is not None:
            self.read_head_config(head_config)
        if calc_config is not None:
            self.read_calc_config(calc_config)

        for _, __ in kwargs.items():
            if _ in self.head_pars:
                self.head_pars[_] = __
            elif _ in self.calc_pars:
                self.calc_pars[_] = __
            else:
                warn(f"Parameter {_} not recognised --- skipped")

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
Вписать delta_T в параметры обработки
"""
