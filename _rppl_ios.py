import json
from warnings import warn

from astropy.io import ascii


class _RPPLSubIOS:
    def __init__(self):
        self.head_pars = None
        self.calc_pars = None
        self.images_table = None
        self.sdfimg_table = None

    def read_head_config(self, filename=None):
        if filename is None:
            filename = f"{self.head_pars['CONFIG_FILENAME']}.json"
        with open(filename, "r") as confile:
            buff_dict = json.load(confile)
            for _ in buff_dict.keys():
                if _ in self.head_pars:
                    self.head_pars[_] = buff_dict[_]
                else:
                    warn(f"Parameter {_} not recognised - skipped")

    def read_calc_config(self, filename=None):
        if filename is None:
            filename = f"{self.calc_pars['CONFIG_FILENAME']}.json"
        with open(filename, "r") as confile:
            buff_dict = json.load(confile)
            for _ in buff_dict.keys():
                if _ in self.calc_pars:
                    self.calc_pars[_] = buff_dict[_]
                else:
                    warn(f"Parameter {_} not recognised - skipped")

    def write_head_config(self, filename=None):
        if filename is None:
            filename = f"{self.head_pars['CONFIG_FILENAME']}.json"
        with open(filename, "w") as confile:
            json.dump(self.head_pars, confile, indent=4)

    def write_calc_config(self, filename=None):
        if filename is None:
            filename = f"{self.calc_pars['CONFIG_FILENAME']}.json"
        with open(filename, "w") as confile:
            json.dump(self.calc_pars, confile, indent=4)

    def read_images_table(self, filename=None):
        if filename is None:
            filename = f"images_table_{self.head_pars['CONFIG_FILENAME']}_" \
                       f"{self.calc_pars['CONFIG_FILENAME']}.txt"
        self.images_table = ascii.read(filename, delimiter="\t",
                                       format="commented_header", fill_values=[(ascii.masked, "nan")])

    def write_images_table(self, filename=None):
        if filename is None:
            filename = f"images_table_{self.head_pars['CONFIG_FILENAME']}_" \
                       f"{self.calc_pars['CONFIG_FILENAME']}.txt"
        ascii.write(self.images_table, filename, overwrite=True, delimiter="\t",
                    format="commented_header", fill_values=[(ascii.masked, "nan")])

    def read_sdfimg_table(self, filename=None):
        if filename is None:
            filename = f"sdfimg_table_{self.head_pars['CONFIG_FILENAME']}_" \
                       f"{self.calc_pars['CONFIG_FILENAME']}.txt"
        self.sdfimg_table = ascii.read(filename, delimiter="\t",
                                       format="commented_header", fill_values=[(ascii.masked, "nan")])

    def write_sdfimg_table(self, filename=None):
        if filename is None:
            filename = f"sdfimg_table_{self.head_pars['CONFIG_FILENAME']}_" \
                       f"{self.calc_pars['CONFIG_FILENAME']}.txt"
        ascii.write(self.sdfimg_table, filename, overwrite=True, delimiter="\t",
                    format="commented_header", fill_values=[(ascii.masked, "nan")])


"""
Стоит внедрить чтение и запись с помощью методов для таблиц
"""
