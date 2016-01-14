import pylab
from desispec.log import get_logger

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

    #   If end is not defined, end may take nfibers value
    #   else it will take max number of fibers value
    if end is None :
        end = f_shape

    if start > end :
        log.error("Warning. Nothing will be plotted since begin value is higher than end value")
   
    #   These values come from a string.
    #   They need this to be treated as integer.
    start   = int(start)
    end     = int(end)
    
    print ivar
    #print len(spectra[0])
    for fiber in xrange(f_shape) :
        if (fiber >= start and fiber <= end) :
            log.info("Plotting fiber %03d" % fiber)
            pylab.plot(spectra[fiber], 'b')
            pylab.plot(ivar[fiber], 'r')
    pylab.show()
