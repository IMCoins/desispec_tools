#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
from desispec_tools.graph_tools         import plot_graph
from desispec_tools.graph_tools         import show_graph
from desispec.log import get_logger

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--frame', type = str, default = None, required = True,
                    help = 'path of frame fits file')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on.\
                    (ex: --fibers=50:60 means that only fibers from 50 to 60\
                    (60 being excluded) will be plotted)')

log         = get_logger()
args        = parser.parse_args()
frame_file  = pyfits.open(args.frame)



# --fibers 1,4,6:8,5

if args.fibers is not None :
    nb = args.fibers.split(',')
    for i in xrange(len(nb)) :
        print 'nb[i]=', nb[i].isdigit()
        if nb[i].isdigit() == False :
            tmp = nb[i].split(':')
            print 'len(tmp)=', len(tmp)
            if ((len(tmp) is 2) and tmp[0].isdigit() == True and tmp[1].isdigit() == True) :
                f_begin = int(tmp[0])
                f_end   = int(tmp[1])
                plot_graph(frame_file, nfibers=None, start=f_begin, end=f_end)
            else :
                log.error("--fibers parsing error.\
                    correct format is either  : --fibers=begin,end (excluded)\
                                      and/or  : --fibers=begin:end (excluded)\
                    You can use : --fibers=2,5,6:8,3,10")
                sys.exit(1)
        else :
            log.info("{DEBUG}")
            nb[i] = int(nb[i])
            plot_graph(frame_file, nfibers=None, start=nb[i], end=None, only=True)
else :
    #   If you did not ask for a specific fiber/group of fiber, they're all plotted.
    plot_graph(frame_file, nfibers=frame_file[0].data.shape[0])


"""

if args.fibers is not None :
    nb = args.fibers.split(':')
    if len(nb) is not 2 or nb[0].isdigit == False or nb[1].isdigit == False  :
        log.error("--fibers parsing error. correct format is : --fibers=begin,end (end is excluded)")
        sys.exit(1)
        # fibers_begin = nb[0]
        # fibers_end   = nb[1]
    else :
        fibers_begin = nb[0]
        fibers_end   = nb[1]
        # fibers_begin = 0
        # fibers_end   = None

"""

# log.info("Showing.......")
# plot_graph(frame_file, nfibers=None, start=fibers_begin, end=fibers_end)


show_graph()
log.info("Script's done")