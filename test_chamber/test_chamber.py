import astropy.io.fits as fits

print("Sandbox code activated")

with fits.open("PaulStar-2021-02-17T14-14-10.fits") as hdu_list:
    hdu_list.verify("fix")
    img_header = hdu_list[0].header
    img_data = hdu_list[0].data

A = img_header["IMAGETYP"]
print(A)

"""

"""
