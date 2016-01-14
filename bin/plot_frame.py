#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
from desispec_tools.graph_tools         import plot_graph
from desispec.log import get_logger

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--frame', type = str, default = None, required = True,
                    help = 'path of frame fits file')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60 means that only fibers from 50 to 60 (60 being excluded) will be plotted)')

log = get_logger()

args = parser.parse_args()


frame_file  = pyfits.open(args.frame)

if args.fibers is not None :
    nb = args.fibers.split(':')
    
    
    if (len(nb) is not 2) or nb[0].isdigit == False or nb[1].isdigit == False  :
         log.error("--fibers parsing error. correct format is : --fibers=begin,end (end is excluded)")
         sys.exit(1)
    fibers_begin = nb[0]
    fibers_end   = nb[1]
else :
    fibers_begin = 0
    fibers_end   = None

log.info("Showing.......")
plot_graph(frame_file, nfibers=None, start=fibers_begin, end=fibers_end)



log.info("Script done")
