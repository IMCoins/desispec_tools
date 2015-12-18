#!/usr/bin/env python

import argparse
from boxcar_extraction import boxcar

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--boot', type = str, default = None, required = True,
                    help = 'path of .fits file to get wavelength from')
parser.add_argument('--pix', type = str, default = None, required = True,
                    help = 'path of .fits file in which specter are interpreted')
args = parser.parse_args()

boxcar(args.boot, args.pix)
