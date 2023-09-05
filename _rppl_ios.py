import json
from warnings import warn

from astropy.io import ascii

from utils import check_out_directory


class _IOHandler:
    def __init__(self):
        self.settings = None
        self.images_table = None

    def read_config(self, filename):
        with open(filename, "r") as config_file:
            buff_dict = json.load(config_file)
            for _ in buff_dict.keys():
                if _ in self.settings:
                    self.settings[_] = buff_dict[_]
                else:
                    warn(f"Parameter {_} not recognised - skipped")
        print(f"Loaded config: {filename}\n")

    def write_config(self, filename):
        with open(filename, "w") as confile:
            json.dump(self.settings, confile, indent=4)

    def read_images_table(self, filename):
        self.images_table = ascii.read(filename, delimiter="\t",
                                       format="commented_header", fill_values=[(ascii.masked, "nan")])

    def write_images_table(self, filename):
        check_out_directory(self.settings["OUT_DIR"])
        ascii.write(self.images_table, self.settings["OUT_DIR"] + filename, overwrite=True, delimiter="\t",
                    format="commented_header", fill_values=[(ascii.masked, "nan")])
        print(f"Saved images table: {filename}\n")
