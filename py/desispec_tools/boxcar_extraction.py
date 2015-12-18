import numpy as np
import astropy.io.fits as pyfits
from numpy.polynomial.legendre import legval, legfit
from desispec.log import get_logger

def u(wave, wavemin, wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

def invert_legendre_polynomial() :
    # input = xmin,xmax,coef_of_x_to_y
    # output = ymin,ymax,coef_of_y_to_x
    pass

################
#   RETURNS FITS FILE INCLUDING ELECTRONS QUANTITY
################




def boxcar(psf, image_file, nfibers=None) :
    """Find and returns  wavelength  spectra and inverse variance

        Parameters
        ----------
        psf     : File Descriptor
        Where the wavelength is collected.

        pix     : File Descriptor
        Interpreted photons using the wavelength.

        

        nfibers : Optional. If left empty, will set the max number of fibers.
                            Max number is 500

        Returns
        spectra
        ivar
        wavelength
        """
    log=get_logger()
    log.info("starting boxcar extraction")
    
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

    # number of pixels in image 
    # we are going to extract one flux per fiber per y pixel (total = nfibers x npix_y)
    npix_y  = flux.shape[0]
    npix_x  = flux.shape[1]

    # wavelength array used to invert 
    # NOT NEEDED WHEN LEGPOL INVERSION DONE IN ROUTINE
    tmp_wave_for_inversion = np.linspace(wavemin, wavemax, 100)
    
    nfibers_to_extract = xcoef.shape[0]
    if nfibers is not None :
        if nfibers>nfibers_to_extract :
            log.warning("only %d fibers will be extracted"%nfibers_to_extract)
        
        nfibers_to_extract = min(nfibers,nfibers_to_extract)
    

    #   Flux as a function of wavelength
    spectra             = np.zeros((nfibers_to_extract,npix_y))
    #   Inverse-variance of spectrum
    spectra_ivar        = np.zeros((nfibers_to_extract,npix_y))
    #   Wavelength
    wave_of_y           = np.zeros((nfibers_to_extract, npix_y))

###
# Using legendre's polynomial to get a spectrum per fiber
###

    for fiber in xrange(nfibers_to_extract) :
        
        log.info("extracting fiber #%03d"%fiber)
        
        #   Determines value of Y, so we can know its coeficient and then its position
        y_of_wave           = legval(u(tmp_wave_for_inversion, wavemin, wavemax), ycoef[fiber])
        coef                = legfit(u(y_of_wave, 0, npix_y),tmp_wave_for_inversion , deg=ycoef[fiber].size)
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
        

    log.info("done boxcar extraction")
    return spectra,spectra_ivar,wave_of_y


