#!/usr/bin/env python

import argparse
import astropy.io.fits as pyfits
from boxcar_extraction import boxcar

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--boot', type = str, default = None, required = True,
                    help = 'path of .fits file to get wavelength from')
parser.add_argument('--pix', type = str, default = None, required = True,
                    help = 'path of .fits file in which specter are interpreted')
args = parser.parse_args()

psf         = pyfits.open(args.boot)
image_file  = pyfits.open(args.pix)

boxcar(psf, image_file)
