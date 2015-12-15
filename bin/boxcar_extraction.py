#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import pylab
from numpy.polynomial.legendre import legval, legfit

def u(wave,wavemin,wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

################
#   RETURNS FITS FILE INCLUDING ELECTRONS QUANTITY
################

def boxcar(file_1, file_2) :
    """Find and returns a fits file thanks to an input spectrum

        Parameters
        ----------
        boot : fits file
        Where the wavelength is collected.
        pix : fits file
        Interpreted photons using the wavelength.
        """
    nfibers = 50
    psf     = pyfits.open(file_1)
#    print psf[0].header
#    print psf.info()
    wavemin = psf[0].header["WAVEMIN"]
    wavemax = psf[0].header["WAVEMAX"]
    xcoef   = psf[0].data
    ycoef   = psf[1].data
    xsig    = psf[2].data       	
#    print xcoef.shape
    image_file  = pyfits.open(file_2)
    flux        = image_file[0].data
    flux_ivar   = image_file[1].data      # inverse de la variance de la valeur des pixels de l'image
    flux_var    = np.zeros(flux_ivar.shape)
    mask        = (flux_ivar > 0)

    flux_var[mask] = 1./flux_ivar[mask]

    npix_y  = flux.shape[0]
    npix_x  = flux.shape[1]
    wave=np.linspace(wavemin,wavemax,10)

    spectra             = np.zeros((nfibers,npix_y)) # Flux as a function of wavelength
    spectra_ivar        = np.zeros((nfibers,npix_y)) # Reverse variance of spectrum
    wave_vs_y_of_fibers = np.zeros((nfibers,npix_y))

# Using legendre's polynomial to get a spectrum per fiber.

    for fiber in xrange(nfibers) :
        print fiber
        y_of_wave   = legval(u(wave, wavemin, wavemax), ycoef[fiber])
        coef        = legfit(u(y_of_wave, 0, npix_y), wave, deg=ycoef[fiber].size)
        wave_of_y   = legval(u(np.arange(npix_y).astype(float), 0, npix_y), coef)
        x_of_y      = legval(u(wave_of_y, wavemin, wavemax), xcoef[fiber])
        x1_of_y     = np.floor(x_of_y).astype(int) - 3
        x2_of_y     = np.floor(x_of_y).astype(int) + 4
        for y in xrange(npix_y) :
            nb_invalidPix=np.sum(flux_ivar[y,x1_of_y[y]:x2_of_y[y]] <= 0) # number of pixels with neg or null flux_ivar
            if nb_invalidPix == 0 :
                spectra[fiber, y] = np.sum(flux[y,x1_of_y[y]:x2_of_y[y]])       # sum of flux
                var               = np.sum(flux_var[y,x1_of_y[y]:x2_of_y[y]])   # sum of variances
                spectra_ivar[fiber, y] = 1./var
        wave_vs_y_of_fibers[fiber]= wave_of_y

    hdulist=pyfits.HDUList([pyfits.PrimaryHDU(spectra),
                            pyfits.ImageHDU(spectra_ivar,name="IVAR"),
                            pyfits.ImageHDU(wave_vs_y_of_fibers,name="WAVELENGTH")])
    # HDU[0]: Spectra
    # HDU[1]: Inverse variance of spectra
    # HDU[2]: Wavelengh
    # HDU[3]: Resolution matrix
    hdulist.writeto("frame.fits",clobber=True)
