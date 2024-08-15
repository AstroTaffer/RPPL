import os
import subprocess
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
from astropy import units as u


def DoAss(file_name, depth=100, sigma=300):
    hdulist = fits.open(file_name, 'update', memmap=False)
    header = hdulist[0].header
    binning = header['XBINNING']
    try:
        cord = SkyCoord(f'{header["OBJRA"]} {header["OBJDEC"]}', unit=(u.hourangle, u.deg), frame='icrs')
    except:
        try:
            cord = SkyCoord(f'{header["ALPHA"]} {header["DELTA"]}', unit=(u.hourangle, u.deg), frame='icrs')
        except:
            print('ERROR: can\'t find coords in header')
            return False
    shell = f'C:/Users/Администратор/AppData/Local/cygwin_ansvr/bin/mintty.exe /bin/bash -l solve-field --ra {cord.ra.deg} --dec {cord.dec.deg} --downsample 2 --radius 0.7 --depth {depth} --sigma {sigma} --no-remove-lines --overwrite --no-verify --no-verify-uniformize --no-verify-dedup -L {0.55 * binning} -H {0.75 * binning} -u app -O -p -r -t 2 -M none -R none -S none -P none -B none -U none -N none -l 60 \"/cygdrive/{file_name.replace(":", "")}\"'
    print(shell)
    try:
        subprocess.check_output(shell, shell=True)
        path = file_name.split('.')
        wcs_path = '.'.join(path[:-1]) + '.wcs'
        axy_path = '.'.join(path[:-1]) + '.axy'
        with open(wcs_path) as f:
            wcs = f.read()
        if os.path.exists(wcs_path):
            os.remove(wcs_path)
        print('Astrometry: solved')
        if os.path.exists(axy_path):
            os.remove(axy_path)
        wcs = wcs.split('HISTORY')[0]
        header_wcs = fits.Header.fromstring(wcs)
        hdulist[0].header = header + header_wcs
    except Exception as e:
        print('Astrometry: error')
        print(e)
        hdulist.close()
        return False
    hdulist.close()
    return True


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