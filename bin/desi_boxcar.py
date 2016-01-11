#!/usr/bin/env python

import argparse
import astropy.io.fits as pyfits
from desispec_tools.boxcar_extraction import boxcar
from desispec_tools.resample import resample_to_same_wavelenght_grid
from desispec_tools.graphic_extraction import show_graph

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf fits file to get wavelength from')
parser.add_argument('-i','--image', type = str, default = None, required = True,
                    help = 'path of image fits file')
parser.add_argument('-o','--outframe', type = str, default = None, required = True,
                    help = 'path of output frame file')
parser.add_argument('-n','--nfibers', type = int, default = None, required = False,
                    help = 'number of fibers (default=all)')
parser.add_argument('--show', action='store_true',
                    help = 'plot result')
parser.add_argument('-r','--resample', action='store_true',
                    help = 'resample to save wavelength grid')


args = parser.parse_args()

psf         = pyfits.open(args.psf)
image_file  = pyfits.open(args.image)

spectra, ivar, wave = boxcar(psf,image_file,args.nfibers)
#boxcar(psf,image_file,args.nfibers)

if args.resample :
    spectra, ivar, wave = resample_to_same_wavelenght_grid(spectra, ivar, wave)


hdulist = pyfits.HDUList([pyfits.PrimaryHDU(spectra),
                        pyfits.ImageHDU(ivar,name="IVAR"),
                        pyfits.ImageHDU(wave,name="WAVELENGTH")])
                        #pyfits.ImageHDU(rdata, name="RESOLUTION")])
hdulist.writeto(args.outframe,clobber=True)

if args.show :
    frame       = pyfits.open(args.outframe)
    show_graph(frame)
