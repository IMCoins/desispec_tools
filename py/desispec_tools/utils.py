import numpy as np
from numpy.polynomial.legendre import legval, legfit


def u(wave, wavemin, wavemax) :
    return 2. * (wave - wavemin)/(wavemax - wavemin) - 1.

def invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, fiber, npix_y, wave_of_y) :
 
    #   Wavelength array used in 'invert_legendre_polynomial'
    wave                = np.linspace(wavemin, wavemax, 100)
    #   Determines value of Y, so we can know its coeficient and then its position
    y_of_wave           = legval(u(wave, wavemin, wavemax), ycoef[fiber])
    coef                = legfit(u(y_of_wave, 0, npix_y), wave, deg=ycoef[fiber].size)
    wave_of_y[fiber]    = legval(u(np.arange(npix_y).astype(float), 0, npix_y), coef)
    #   Determines wavelength intensity (x) based on Y
    x_of_y              = legval(u(wave_of_y[fiber], wavemin, wavemax), xcoef[fiber])
    return x_of_y
