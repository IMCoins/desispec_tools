#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import pylab
from numpy.polynomial.legendre import legval, legfit

def u(wave, wavemin, wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

################
#   RETURNS FITS FILE INCLUDING ELECTRONS QUANTITY
################

def boxcar(psf, image_file, graph=0) :
    """Find and returns a fits file thanks to an input wavelength and spectrum

        Parameters
        ----------
        psf     : File Descriptor
        Where the wavelength is collected.

        pix     : File Descriptor
        Interpreted photons using the wavelength.

        graph   : Optional. If empty, will return a fits file.
                            If not and set to an int > 0, will show a graph
                            of the spectrum.
        """
    nfibers = 500
    wavemin = psf[0].header["WAVEMIN"]
    wavemax = psf[0].header["WAVEMAX"]
    xcoef   = psf[0].data
    ycoef   = psf[1].data
    xsig    = psf[2].data
    flux        = image_file[0].data
    flux_ivar   = image_file[1].data
    flux_var    = np.zeros(flux_ivar.shape)

    mask        = (flux_ivar > 0)
    flux_var[mask] = 1./flux_ivar[mask]

    npix_y  = flux.shape[0]
    npix_x  = flux.shape[1]
    wave    = np.linspace(wavemin, wavemax, 10)

    spectra             = np.zeros((nfibers,npix_y))
    spectra_ivar        = np.zeros((nfibers,npix_y))
    wave_of_y           = np.zeros((nfibers, npix_y))

###
# Using legendre's polynomial to get a spectrum per fiber
###

    for fiber in xrange(nfibers) :
        print fiber
        y_of_wave           = legval(u(wave, wavemin, wavemax), ycoef[fiber])
        coef                = legfit(u(y_of_wave, 0, npix_y), wave, deg=ycoef[fiber].size)
        wave_of_y[fiber]    = legval(u(np.arange(npix_y).astype(float), 0, npix_y), coef)
        x_of_y              = legval(u(wave_of_y[fiber], wavemin, wavemax), xcoef[fiber])
        x1_of_y             = np.floor(x_of_y).astype(int) - 3
        x2_of_y             = np.floor(x_of_y).astype(int) + 4
        for y in xrange(npix_y) :
            nb_invalidPix   = np.sum(flux_ivar[y, x1_of_y[y]:x2_of_y[y]] <= 0)
            if nb_invalidPix == 0 :
                spectra[fiber, y]   = np.sum(flux[y, x1_of_y[y]:x2_of_y[y]])
                var                 = np.sum(flux_var[y, x1_of_y[y]:x2_of_y[y]])
                spectra_ivar[fiber, y] = 1./var
        if graph :
            pylab.plot(spectra[fiber])

    if graph == 0 :
        hdulist=pyfits.HDUList([pyfits.PrimaryHDU(spectra),
                                pyfits.ImageHDU(spectra_ivar,name="IVAR"),
                                pyfits.ImageHDU(wave_of_y,name="WAVELENGTH")])
                                #pyfits.ImageHDU(rdata, name="RESOLUTION")])
        hdulist.writeto("frame.fits",clobber=True)
        # hdulist[0]: Spectra
        # hdulist[1]: Inverse variance of spectra
        # hdulist[2]: Wavelengh
        # hdulist[3]: Resolution matrix
    else :
        pylab.show()
