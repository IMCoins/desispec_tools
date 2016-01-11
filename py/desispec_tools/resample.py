import numpy as np

from desispec.interpolation import resample_flux

def resample_to_same_wavelenght_grid(spectra, ivar, wave) :

    # chose the average wavelength of all fibers
    same_wave=np.mean(wave,axis=0)

    nfibers=spectra.shape[0]

    # declare output
    resampled_spectra = np.zeros((nfibers,same_wave.size))
    resampled_ivar = np.zeros((nfibers,same_wave.size))

    # loop on fibers
    for fiber in range(nfibers) :        
        resampled_spectra[fiber],resampled_ivar[fiber] = resample_flux(same_wave,wave[fiber],spectra[fiber],ivar[fiber])
    
    return resampled_spectra, resampled_ivar, same_wave


