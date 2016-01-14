import pylab
from desispec.log import get_logger
import sys
import numpy as np

def plot_graph(frame, nfibers=None, start=0, end=None) :
    """Plot graph from a given spectra from a fits file.

    ----------
    Parameters
    ----------

    frame : File Directory
    Where the spectra is collected to be plot.

    nfibers : Optional. If left empty, will set the max number of fibers.

    start   : Optional. If left empty, will be set to 0.
    end     : Optional. If left empty, will be set to nfibers.
    Define start/end argument to choose from frame parameter which specific fiber to plot.
    Example     -> plot_graph(frame, 500, 10, 20) will take 500 fibers into account, and plot from the 10th fiber to the 20th (included).
    """

    log         = get_logger()
    spectra     = frame[0].data
    ivar        = frame[1].data
    wave        = frame[2].data
    f_shape     = spectra.shape[0]
    
    if nfibers is not None :
        if nfibers > f_shape :
            log.warning("only %d fibers will be shown" % f_shape)
        f_shape = min(nfibers, f_shape)
    
    #   These values come from a string.
    #   They need this to be treated as integer.
    start   = int(start)
    end     = int(end)

    #   If end is not defined, end may take nfibers value
    #   else it will take max number of fibers value
    if end is None :
        end = f_shape
    elif end > f_shape :
        log.warning("will plot only %d fibers"%f_shape)
        end = f_shape
    
    if start > end :
        log.error("Warning. Nothing will be plotted since begin value is higher than end value")
   
   
    
    show_errors = True
    
    if show_errors :
        err = np.sqrt(1./(ivar+(ivar==0)))*(ivar>0)
        
    
    #print len(spectra[0])
    for fiber in xrange(start,end) :
        log.info("Plotting fiber %03d" % fiber)
            
        if show_errors :
                
            if len(wave.shape)>1 :
                pylab.errorbar(wave[fiber],spectra[fiber],err[fiber],fmt="o-")
            else :
                pylab.errorbar(wave,spectra[fiber],err[fiber],fmt="o-")
        else :
            if len(wave.shape)>1 :
                pylab.plot(wave[fiber],spectra[fiber],"-")
            else :
                pylab.plot(wave,spectra[fiber],"-")
    
    pylab.xlabel("Wavelength [A]")


    show_2d=True

    if show_2d :
        pylab.figure()
        if len(wave.shape)==1 :
            print wave[0],wave[-1],start-0.5,end-0.5
            print spectra[start:end].shape
            pylab.imshow(spectra[start:end],aspect=spectra.shape[1]/float(end-start),extent=(wave[0],wave[-1],start-0.5,end-0.5),origin=0.,interpolation="nearest")
            pylab.xlabel("Wavelength [A]")
            pylab.ylabel("Fiber #")
        else :
            pylab.imshow(spectra[start:end],aspect=spectra.shape[1]/float(end-start),extent=(0,spectra.shape[1],start-0.5,end-0.5),origin=0.,interpolation="nearest")
            pylab.xlabel("Y CCD")
            pylab.ylabel("Fiber #")
        pylab.colorbar()
            


    pylab.show()
