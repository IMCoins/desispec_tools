import pylab

def show_graph(frame, nfibers=None) :
    """Shows graph from a given spectra from a fits file.

    ----------
    Parameters
    ----------

    frame : File Directory
    Where the spectra is collected

    nfibers : Optional. If left empty, will set the max number of fibers.
    """

    spectra     = frame[0].data
    f_shape     = spectra.shape[0]

    if nfibers is not None :
        if nfibers > f_shape :
#            log.warning("only %d fibers will be shown" % f_shape)
            print 'WARNING'
        nfibers = min(nfibers, f_shape)

    for fiber in xrange(f_shape) :
        pylab.plot(spectra[fiber])
    pylab.show()
