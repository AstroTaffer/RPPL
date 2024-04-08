import glob
import os
import subprocess

import astropy.io.fits as fits
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std
from astroquery.astrometry_net import AstrometryNet
from photutils.detection import DAOStarFinder
from astropy import units as u


def Astrometry(path, files, C):
    ast = AstrometryNet()
    ast.api_key = 'hipfhzhlzygnlvix'
    path2save = path + r'\AstrometryDone'
    if not os.path.exists(path2save):
        os.mkdir(path2save)
    for count, item in enumerate(files):
        try:
            name = item.split('\\')[-1]
            print(f"Astrometry ({count}/{len(files)}): working on {name}")
            # read file, copy data and header
            with fits.open(item, mode='update', memmap=False) as hdulist:
                # hdulist = fits.open(item, 'update', memmap=False)
                # del hdulist[0].header['COMMENT']
                Header = hdulist[0].header
                try:
                    buf = Header['CD1_1']
                    hdulist.close()
                    # print('done')
                    continue
                except:
                    pass
                Data = hdulist[0].data.copy()
                hdulist.verify('fix')
                hdulist.close()
                # gaussian convolution
                kernel = Gaussian2DKernel(x_stddev=1)
                Data = convolve(Data, kernel)
                # extract background
                Data -= np.median(Data)
                Bkg_sigma = mad_std(Data)
                # # mask bad row
                # mask = np.zeros(Data.shape, dtype=bool)
                # mask[90:110, 0:Data.shape[1]] = True
                daofind = DAOStarFinder(fwhm=4.5, threshold=5. * Bkg_sigma, sharplo=0.25)
                # Sources = daofind(Data, mask=mask)
                Sources = daofind(Data)
                # print(Sources.info)
                # plt.imshow(Data, cmap=cm.Greys_r, aspect='equal',
                #            norm=Normalize(vmin=-30, vmax=150), interpolation='nearest')
                # plt.scatter(Sources['xcentroid'], Sources['ycentroid'], s=40, facecolors='none', edgecolors='r')
                # plt.show()
                # Sort sources in ascending order
                Sources.sort('flux')
                Sources.reverse()
                # ast.show_allowed_settings()
                image_width = Header['NAXIS2']
                image_height = Header['NAXIS1']
                # print(Sources)
                wcs_header = ast.solve_from_source_list(Sources['xcentroid'],
                                                        Sources['ycentroid'],
                                                        image_width, image_height,
                                                        solve_timeout=120,
                                                        center_ra=C.ra.degree,
                                                        center_dec=C.dec.degree,
                                                        radius=0.25,
                                                        downsample_factor=2,
                                                        scale_lower=1.2,
                                                        scale_upper=1.4,
                                                        scale_units='arcsecperpix'
                                                        )
                # hdu = fits.PrimaryHDU(hdulist[0].data, Header + wcs_header)
                hdulist[0].header = Header + wcs_header
                # hdulist.close()
                # hdu.writeto(path2save + '\\' + name, overwrite=True)
                hdulist.writeto(item, overwrite=True)
                print('done')
        except Exception as e:
            print(e)
    return path


class DoAstrometry:

    def __init__(self, key):
        self.ast = AstrometryNet()
        self.ast.api_key = key

    def compute(self, input_path):
        try:
            # read file, copy data and header
            hdulist = fits.open(input_path, 'update', memmap=False)
            Header = hdulist[0].header
            try:
                buf = Header['CD1_1']
                hdulist.close()
                return True
            except:
                pass
            C = SkyCoord(f'{Header["OBJRA"]} {Header["OBJDEC"]}', unit=(u.hourangle, u.deg), frame='icrs')
            Data = hdulist[0].data.copy()
            hdulist.verify('fix')
            hdulist.close()
            # gaussian convolution
            kernel = Gaussian2DKernel(x_stddev=1)
            Data = convolve(Data, kernel)
            # extract background
            Data -= np.median(Data)
            Bkg_sigma = mad_std(Data)
            # # mask bad row
            # mask = np.zeros(Data.shape, dtype=bool)
            # mask[90:110, 0:Data.shape[1]] = True
            daofind = DAOStarFinder(fwhm=4.5, threshold=5. * Bkg_sigma, sharplo=0.25)
            # Sources = daofind(Data, mask=mask)
            Sources = daofind(Data)
            # Sort sources in ascending order
            Sources.sort('flux')
            Sources.reverse()
            # ast.show_allowed_settings()
            image_width = Header['NAXIS2']
            image_height = Header['NAXIS1']
            # print(Sources)
            wcs_header = self.ast.solve_from_source_list(Sources['xcentroid'],
                                                         Sources['ycentroid'],
                                                         image_width, image_height,
                                                         solve_timeout=30,
                                                         center_ra=C.ra.degree,
                                                         center_dec=C.dec.degree,
                                                         radius=0.25,
                                                         downsample_factor=2,
                                                         scale_lower=1.2,
                                                         scale_upper=1.4,
                                                         scale_units='arcsecperpix'
                                                         )
            hdulist[0].header = Header + wcs_header
            hdulist.writeto(input_path, overwrite=True)
            hdulist.close()
            return True
        except Exception as e:
            print(e)
            return False


    # WCS_success, WCS_Name = Astrometry('/dev/shm/Proc_Temp.fits', \
    #                                    Center_Ra, Center_Dec, 100, 300, 2, Header['XBINNING'])
def DoAss(file_name, Ra, Dec, Depth, Sigma, sip, binning):
    # success = False
    New_name = file_name.replace('.fits', '_WCS.fits')

    shell = r'C:\cygwin64\Cygwin.bat solve-field'
    shell = shell + ' --ra '
    shell = shell + str(Ra)
    shell = shell + ' --dec '
    shell = shell + str(Dec)
    shell = shell + ' --downsample 2'
    shell = shell + ' --radius 0.5'
    shell = shell + ' --depth ' + str(Depth)
    shell = shell + ' --sigma ' + str(Sigma)
    shell = shell + ' --no-remove-lines --no-verify-uniformize --no-verify-dedup'
    shell = shell + ' -L '
    shell = shell + str(0.55*binning)
    shell = shell + ' -H '
    shell = shell + str(0.75*binning)
    shell = shell + ' -u app -O -p -r '
    if sip == 0:
        shell = shell + '-T'
    else:
        shell = shell + '-t ' + str(int(sip))
    shell = shell + ' -W none --no-verify --crpix-center -M none -R none -O -l 60'  #--guess-scale  --no-verify
    shell = shell + ' -S none -B none -U xylist-indx.xyls -N '  # -i none
    # shell = shell + ' -S none -B none --temp-axy -U xylist-indx.xyls -N none'  # -i none
    shell = shell + '\"' + New_name + '\"'
    shell = shell + ' ' + '\"' + file_name + '\"' + '\n'
    print(shell)
    devnull = open(os.devnull, 'w')
    # returncode = subprocess.call(shell, stdout=devnull, stderr=devnull, shell=True)
    returncode = subprocess.check_output(shell, shell=True, stderr=devnull)

    if returncode == 0 and os.path.isfile(New_name):
        #         print('Astrometry: solved')
        return True
    else:
        # clean tmp files
        tmp_files = glob.glob('/tmp/tmp.*')
        for f in tmp_files:
            os.remove(f)
    #         print ('Astrometry: error')

    #     except:
    #         print ('Astrometry: error')
    #         pass
    return False


"""
Options include:
  -h / --help: print this help message
  -v / --verbose: be more chatty -- repeat for even more verboseness
  -D / --dir <directory>: place all output files in the specified directory
  -o / --out <base-filename>: name the output files with this base name
  -b / --backend-config <filename>: use this config file for the "backend"
          program
  --backend-batch: run backend once, rather than once per input file
  -f / --files-on-stdin: read filenames to solve on stdin, one per line
  -p / --no-plots: don't create any plots of the results
  --plot-scale <scale>: scale the plots by this factor (eg, 0.25)
  --plot-bg <filename (JPEG)>: set the background image to use for plots
  -G / --use-wget: use wget instead of curl
  -O / --overwrite: overwrite output files if they already exist
  -K / --continue: don't overwrite output files if they already exist; continue
          a previous run
  -J / --skip-solved: skip input files for which the 'solved' output file
          already exists; NOTE: this assumes single-field input files
  -N / --new-fits <filename>: output filename of the new FITS file containing
          the WCS header; "none" to not create this file
  -Z / --kmz <filename>: create KMZ file for Google Sky.  (requires wcs2kml)
  -i / --scamp <filename>: create image object catalog for SCAMP
  -n / --scamp-config <filename>: create SCAMP config file snippet
  -U / --index-xyls <filename>: output filename for xylist containing the image
          coordinate of stars from the index
  --just-augment: just write the augmented xylist files; don't run backend.
  -7 / --no-delete-temp: don't delete temp files (for debugging)

  -L / --scale-low <scale>: lower bound of image scale estimate
  -H / --scale-high <scale>: upper bound of image scale estimate
  -u / --scale-units <units>: in what units are the lower and upper bounds?
     choices:  "degwidth", "degw", "dw"   : width of the image, in degrees (default)
               "arcminwidth", "amw", "aw" : width of the image, in arcminutes
               "arcsecperpix", "app": arcseconds per pixel
  -8 / --parity <pos/neg>: only check for matches with positive/negative parity
          (default: try both)
  -c / --code-tolerance <distance>: matching distance for quads (default: 0.01)
  -E / --pixel-error <pixels>: for verification, size of pixel positional error
          (default: 1)
  -q / --quad-size-min <fraction>: minimum size of quads to try, as a fraction
          of the smaller image dimension, default: 0.1
  -Q / --quad-size-max <fraction>: maximum size of quads to try, as a fraction
          of the image hypotenuse, default 1.0
  --odds-to-tune-up <odds>: odds ratio at which to try tuning up a match that
          isn't good enough to solve (default: 1e6)
  --odds-to-solve <odds>: odds ratio at which to consider a field solved
          (default: 1e9)
  --odds-to-reject <odds>: odds ratio at which to reject a hypothesis (default:
          1e-100)
  --odds-to-stop-looking <odds>: odds ratio at which to stop adding stars when
          evaluating a hypothesis (default: HUGE_VAL)
  --use-sextractor: use SExtractor rather than built-in image2xy to find sources
  --sextractor-config <filename>: use the given SExtractor config file (default:
          etc/sextractor.conf).  Note that CATALOG_NAME and CATALOG_TYPE values
          will be over-ridden by command-line values.  This option implies
          --use-sextractor.
  --sextractor-path <filename>: use the given path to the SExtractor executable.
          Default: just 'sex', assumed to be in your PATH.  Note that you can
          give command-line args here too (but put them in quotes), eg:
          --sextractor-path 'sex -DETECT_TYPE CCD'.  This option implies
          --use-sextractor.
  -3 / --ra <degrees or hh:mm:ss>: only search in indexes within 'radius' of the
          field center given by 'ra' and 'dec'
  -4 / --dec <degrees or [+-]dd:mm:ss>: only search in indexes within 'radius'
          of the field center given by 'ra' and 'dec'
  -5 / --radius <degrees>: only search in indexes within 'radius' of the field
          center given by ('ra', 'dec')
  -d / --depth <number or range>: number of field objects to look at, or range
          of numbers; 1 is the brightest star, so "-d 10" or "-d 1-10" mean look
          at the top ten brightest stars only.
  --objs <int>: cut the source list to have this many items (after sorting, if
          applicable).
  -l / --cpulimit <seconds>: give up solving after the specified number of
          seconds of CPU time
  -r / --resort: sort the star brightnesses by background-subtracted flux; the
          default is to sort using acompromise between background-subtracted and
          non-background-subtracted flux
  -6 / --extension <int>: FITS extension to read image from.
  -2 / --no-fits2fits: don't sanitize FITS files; assume they're already valid
  --invert: invert the image (for black-on-white images)
  -z / --downsample <int>: downsample the image by factor <int> before running
          source extraction
  --no-background-subtraction: don't try to estimate a smoothly-varying sky
          background during source extraction.
  --sigma <float>: set the noise level in the image
  -9 / --no-remove-lines: don't remove horizontal and vertical overdensities of
          sources.
  --uniformize <int>: select sources uniformly using roughly this many boxes
          (0=disable; default 10)
  --no-verify-uniformize: don't uniformize the field stars during verification
  --no-verify-dedup: don't deduplicate the field stars during verification
  -0 / --no-fix-sdss: don't try to fix SDSS idR files.
  -C / --cancel <filename>: filename whose creation signals the process to stop
  -S / --solved <filename>: output file to mark that the solver succeeded
  -I / --solved-in <filename>: input filename for solved file
  -M / --match <filename>: output filename for match file
  -R / --rdls <filename>: output filename for RDLS file
  --sort-rdls <column>: sort the RDLS file by this column; default is ascending;
          use "-column" to sort "column" in descending order instead.
  --tag <column>: grab tag-along column from index into RDLS file
  --tag-all: grab all tag-along columns from index into RDLS file
  -j / --scamp-ref <filename>: output filename for SCAMP reference catalog
  -B / --corr <filename>: output filename for correspondences
  -W / --wcs <filename>: output filename for WCS file
  -P / --pnm <filename>: save the PNM file as <filename>
  -k / --keep-xylist <filename>: save the (unaugmented) xylist to <filename>
  -A / --dont-augment: quit after writing the unaugmented xylist
  -V / --verify <filename>: try to verify an existing WCS file
  -y / --no-verify: ignore existing WCS headers in FITS input images
  -g / --guess-scale: try to guess the image scale from the FITS headers
  --crpix-center: set the WCS reference point to the image center
  --crpix-x <pix>: set the WCS reference point to the given position
  --crpix-y <pix>: set the WCS reference point to the given position
  -T / --no-tweak: don't fine-tune WCS by computing a SIP polynomial
  -t / --tweak-order <int>: polynomial order of SIP WCS corrections
  -m / --temp-dir <dir>: where to put temp files, default /tmp
The following options are valid for xylist inputs only:
  -F / --fields <number or range>: the FITS extension(s) to solve, inclusive
  -w / --width <pixels>: specify the field width
  -e / --height <pixels>: specify the field height
  -X / --x-column <column-name>: the FITS column containing the X coordinate of
          the sources
  -Y / --y-column <column-name>: the FITS column containing the Y coordinate of
          the sources
  -s / --sort-column <column-name>: the FITS column that should be used to sort
          the sources
  -a / --sort-ascending: sort in ascending order (smallest first); default is
          descending order

Note that most output files can be disabled by setting the filename to "none".
 (If you have a sick sense of humour and you really want to name your output
  file "none", you can use "./none" instead.)
"""