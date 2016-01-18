import numpy as np
from desispec.log import get_logger



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


	flux        = image_file[0].data
	flux_ivar   = image_file[1].data

    #	Usual nfibers' check.
    nfibers_to_extract = xcoef.shape[0]
    if nfibers is not None :
        if nfibers > nfibers_to_extract :
           log.warning("only %d fibers will be extracted" % nfibers_to_extract)
        nfibers_to_extract = min(nfibers,nfibers_to_extract)

    wave_of_y           = np.zeros((nfibers_to_extract, npix_y))
    
    #	Pixel's number in image.
	npix_y  = flux.shape[0]
	for fiber in xrange(nfibers_to_extract) :
		x1_of_y, x2_of_y = invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, fiber, npix_y, wave_of_y)
		for y in xrange(npix_y) :
			#	Checking if there's a dead pixel
            nb_invalidPix   = np.sum(flux_ivar[y, x1_of_y[y]:x2_of_y[y]] <= 0)
            if nb_invalidPix == 0 :
				x = np.arange(x1_of_y,x2_of_y)
				prof = np.exp(-(x - xc) ** 2 / 2 / sigmax ** 2)
				prof /= np.sum(prof)
				spectrum_ivar		= np.sum(flux_ivar[y,x1:x2] * prof ** 2)
				spectrum[fiber, y] 	= np.sum(flux_ivar[y,x1:x2] * image[y,x1:x2] * prof) / spectrum_ivar

	log.info("Profile extraction's done")