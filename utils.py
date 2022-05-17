import astropy.io.fits as fits


class RPPLUtils:
    @staticmethod
    def read_fits_file(file_path):
        with fits.open(file_path) as hdu_list:
            hdu_list.verify("fix")
            img_header = hdu_list[0].header
            img_data = hdu_list[0].data

        return img_header, img_data
