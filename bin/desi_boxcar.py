#!/usr/bin/env python

import argparse
import astropy.io.fits as pyfits
from desispec_tools.boxcar_extraction import boxcar

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf fits file to get wavelength from')
parser.add_argument('-i','--image', type = str, default = None, required = True,
                    help = 'path of image fits file')
parser.add_argument('-o','--outframe', type = str, default = None, required = True,
                    help = 'path of output frame file')
parser.add_argument('-n','--nfibers', type = int, default = None, required = False,
                    help = 'number of fibers (default=all)')


args = parser.parse_args()

psf         = pyfits.open(args.psf)
image_file  = pyfits.open(args.image)

spectra,ivar,wave = boxcar(psf,image_file,args.nfibers)

hdulist=pyfits.HDUList([pyfits.PrimaryHDU(spectra),
                        pyfits.ImageHDU(ivar,name="IVAR"),
                        pyfits.ImageHDU(wave,name="WAVELENGTH")])
                                #pyfits.ImageHDU(rdata, name="RESOLUTION")])
hdulist.writeto(args.outframe,clobber=True)

