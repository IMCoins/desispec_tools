import numpy as np
import astropy.io.fits as pyfits
import pylab
from numpy.polynomial.legendre import legval, legfit

def u(wave, wavemin, wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

################
#   RETURNS FITS FILE INCLUDING ELECTRONS QUANTITY
################

def boxcar(psf, image_file, graph=0, nfibers=500) :
    """Find and returns a fits file thanks to an input wavelength and spectrum

        Parameters
        ----------
        psf     : File Descriptor
        Where the wavelength is collected.

        pix     : File Descriptor
        Interpreted photons using the wavelength.

        graph   : Optional. If left empty, will return a fits file.
                            If not and set to an int > 0, will show a graph
                            of the spectrum.

        nfibers : Optional. If left empty, will set the max number of fibers.
                            Max number is 500
        """
    wavemin = psf[0].header["WAVEMIN"]
    wavemax = psf[0].header["WAVEMAX"]
    xcoef   = psf[0].data
    ycoef   = psf[1].data
    xsig    = psf[2].data

    flux        = image_file[0].data
    #   Inverse variance of the image's value
    flux_ivar   = image_file[1].data
    #   Variance based on inverse variance's size
    flux_var    = np.zeros(flux_ivar.shape)

    #   Applying a mask that keeps positive value to get the Variance by inversing the inverse variance.
    mask        = (flux_ivar > 0)
    flux_var[mask] = 1./flux_ivar[mask]

    #
    npix_y  = flux.shape[0]
    npix_x  = flux.shape[1]
    wave    = np.linspace(wavemin, wavemax, 10)

    #   Flux as a function of wavelength
    spectra             = np.zeros((nfibers,npix_y))
    #   Inverse-variance of spectrum
    spectra_ivar        = np.zeros((nfibers,npix_y))
    #   Wavelength
    wave_of_y           = np.zeros((nfibers, npix_y))

###
# Using legendre's polynomial to get a spectrum per fiber
###

    for fiber in xrange(nfibers) :
        print fiber
        #   Determines value of Y, so we can know its coeficient and then its position
        y_of_wave           = legval(u(wave, wavemin, wavemax), ycoef[fiber])
        coef                = legfit(u(y_of_wave, 0, npix_y), wave, deg=ycoef[fiber].size)
        wave_of_y[fiber]    = legval(u(np.arange(npix_y).astype(float), 0, npix_y), coef)
        #   Determines wavelength intensity (x) based on Y
        x_of_y              = legval(u(wave_of_y[fiber], wavemin, wavemax), xcoef[fiber])
        #   Ascertain X by using low and high uncertainty
        x1_of_y             = np.floor(x_of_y).astype(int) - 3
        x2_of_y             = np.floor(x_of_y).astype(int) + 4
        for y in xrange(npix_y) :
            #   Checking if there's a dead pixel
            nb_invalidPix   = np.sum(flux_ivar[y, x1_of_y[y]:x2_of_y[y]] <= 0)
            if nb_invalidPix == 0 :
                #   Sum of flux
                spectra[fiber, y]   = np.sum(flux[y, x1_of_y[y]:x2_of_y[y]])
                #   Sum of variance
                var                 = np.sum(flux_var[y, x1_of_y[y]:x2_of_y[y]])
                #   Spectrum of inverse variance
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
