from warnings import warn

from _rppl_ios import _IOHandler
from _rppl_ffl import _Lister
from _rppl_clb import _Calibrator


class RPPLData(_IOHandler, _Lister, _Calibrator):
    # noinspection PyMissingConstructor
    def __init__(self, config="default_config.json", **kwargs):
        self.images_table = None
        self.settings = {
            "IN_DIRS": None,
            "CALIBRATION_TYPE": None,  # DO_DARK | DO_FLAT | DO_BOTH
            "OUT_DIR": None,
            "EPSILON_TEMP": None,  # [C]
            "MIN_IMAGES_FOR_MASTER": None,
            "CALIBRATED_BITPIX": None}  # 16 | -64

        # Get config settings
        self.read_config(config)

        # Get kwargs settings
        for _, __ in kwargs.items():
            if _ in self.settings:
                self.settings[_] = __
            else:
                warn(f"Parameter {_} not recognised - skipped")

        # Check for unset settings
        for _ in self.settings.keys():
            if self.settings[_] is None:
                warn(f"Parameter {_} not set")

        print(self.settings)
