import numpy as np
from desispec.log import get_logger
from desispec_tools.utils import invert_legendre_polynomial



def profile_extraction(psf, image_file, nfibers=None) :
	"""
	"""
	log=get_logger()
        log.info("Starting profile extraction...")
    
        wavemin = psf[0].header["WAVEMIN"]
        wavemax = psf[0].header["WAVEMAX"]
        xcoef   = psf[0].data
        ycoef   = psf[1].data
        xsig    = psf[2].data
        print "xsig.shape=",xsig.shape
        print "xsig=",xsig
        
	flux      = image_file[0].data
	flux_ivar = image_file[1].data

        #	Usual nfibers' check.
        nfibers_to_extract = xcoef.shape[0]
        if nfibers is not None :
                if nfibers > nfibers_to_extract :
                        log.warning("only %d fibers will be extracted" % nfibers_to_extract)
                nfibers_to_extract = min(nfibers,nfibers_to_extract)

        #	Pixel's number in image.
        npix_y  = flux.shape[0]
        
        wave_of_y = np.zeros((nfibers_to_extract, npix_y))
        spectrum  = np.zeros((nfibers_to_extract, npix_y))
        spectrum_ivar = np.zeros((nfibers_to_extract, npix_y))
        
	npix_y  = flux.shape[0]
	for fiber in xrange(nfibers_to_extract) :
                # inversion to get x coordinate of trace given y
                xc_of_y = invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, fiber, npix_y, wave_of_y)
                # range of x
                x1_of_y             = np.floor(xc_of_y).astype(int) - 3
                x2_of_y             = np.floor(xc_of_y).astype(int) + 4
		for y in xrange(npix_y) :
                        x1 = x1_of_y[y]
                        x2 = x2_of_y[y]
                        xc = xc_of_y[y]
                        x = np.arange(x1,x2) 
                        # slow and incorrect because we want the integral of the psf in the pixel
                        prof = np.exp(-(x - xc) ** 2 / 2 / xsig[fiber] ** 2)
                        prof /= np.sum(prof)
                        spectrum_ivar[fiber, y]		= np.sum(flux_ivar[y,x1:x2] * prof ** 2)
                        if spectrum_ivar[fiber, y] > 0 :
                                spectrum[fiber, y] 	= np.sum(flux_ivar[y,x1:x2] * flux[y,x1:x2] * prof) / spectrum_ivar[fiber, y]

	log.info("Profile extraction's done")
        return spectrum,spectrum_ivar,wave_of_y
