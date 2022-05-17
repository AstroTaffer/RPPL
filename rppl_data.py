from warnings import warn
from os import listdir

from _rppl_ios import _RPPLSubIOS
from utils import RPPLUtils


class RPPLData(_RPPLSubIOS):
    # noinspection PyMissingConstructor
    def __init__(self, **kwargs):
        self.pars = {"images_directory": None,
                     "images_filter": None,
                     "config_file": kwargs.get("config_file"),

                     "autosave_calcs": True,
                     "autosave_plots": True}

        if self.pars["config_file"] is not None:
            self.read_config()

        for _, __ in kwargs.items():
            if _ in self.pars:
                self.pars[_] = __
            else:
                print(f"Parameter {_} not recognised --- skipped")

        for _ in self.pars.keys():
            if self.pars[_] is None:
                warn(f"Parameter {_} not set")

        if self.pars["autosave_calcs"]:
            self.write_config()

        self.images_namelist = {"Object": [],
                                "Bias": [],
                                "Dark": [],
                                "Flat": []}

        dir_content = listdir(self.pars["images_directory"])
        _ = 0
        while _ < len(dir_content):
            image_name = dir_content[_]
            if not (image_name.count(".fits") or image_name.count(".fit") or image_name.count(".fts")):
                dir_content.remove(image_name)
                continue
            # noinspection PyTypeChecker
            image_header = RPPLUtils.read_fits_file(self.pars["images_directory"] + image_name)[0]
            if not image_header["FILTER"] == self.pars["images_filter"]:
                dir_content.remove(image_name)
            else:
                self.images_namelist[image_header["IMAGETYP"]].append(image_name)
            _ += 1
        print(f"{len(dir_content)} images in total")
        print(f"{len(self.images_namelist['Object'])} Object")
        print(f"{len(self.images_namelist['Bias'])} Bias")
        print(f"{len(self.images_namelist['Dark'])} Dark")
        print(f"{len(self.images_namelist['Flat'])} Flat")


"""

"""
