#!/usr/bin/env python

################################################################################
#                                                                              # 
#  PLOTS SOME QUANTITIES CALCULATED IN-SITU                                    # 
#                                                                              # 
################################################################################

from __future__ import print_function
import sys; sys.dont_write_bytecode = True
import numpy as np
import matplotlib as mpl
mpl.use('agg')
from  matplotlib import pyplot as plt
mpl.rcParams.update({'font.size':18})

def plot_q_of_t(qname   = 'mdot_eh',
                diag    = None,
                path    = None,
                op      = None,
                save    = False,
                display = False,
                ymin    = None,
                ymax    = None,
                log     = False):
    """Plot quantity q, given in the diag dict object by qname,
    as a function of time. Options are:

    qname -- String, Name of quantity to plot out of diag. Default is 'modt_eh'

    diag  -- String, diag object attained by using hdf5_to_dict.load_diag.
             Default is None, but must be set unless path is not None.
          
    path  -- String, path to directory containing diag file. Default is None,
             but must be set unless diag is not None.
          
    op    -- unary function, operation to apply to q before plotting.
             default is the identity.
          
    save  -- Save the plot to a file? Plot is saved into
            current working directory. Default is False.

    ymin -- Minimum y value. Defaults to letting the code choose.

    ymax -- Maximum y value. Defaults to letting the code choose.

    display -- Whether or not to show the plot. Default is False.    

    """
    if diag is None:
        if path is None:
            raise ValueError("Either path or diag must be non-none.")
        import hdf5_to_dict as io
        diag = io.load_diag(path)
    t = diag['t']
    q = diag[qname]
    if op is not None:
        q = op(q)
    if log:
        plt.loglog(t,q)
    else:
        plt.plot(t,q)
    plt.xlabel('t (code units)')
    plt.ylabel(qname + ' (code units)')
    if ymin is not None or ymax is not None:
        if ymin is None:
            ymin = np.min(q)
        if ymax is None:
            ymin = np.max(q)
        plt.ylim(ymin,ymax)
    if save:
        for postfix in ['pdf', 'png']:
            name = '{}_{}.{}'.format(qname,
                                     'log' if log else 'lin',
                                     postfix)
            plt.savefig(name, bbox_inches='tight')
    if display:
        plt.show()
    return

if __name__ =="__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot an in-situ quantity in bhlight.')
    parser.add_argument('dumpdir',
                        type=str,
                        help='Path to directory containing dumps and diag file.')
    parser.add_argument('-q','--qname',
                        dest='qname',
                        type=str,
                        help='Name of quantity to plot.',
                        default='mdot_eh')
    parser.add_argument('-i','--invert',
                        action='store_true',
                        dest='invert',
                        help='Invert the y-axis of the plot.')
    parser.add_argument('-l','--log',
                        action='store_true',
                        dest='log',
                        help='Make y axis on a log scale')
    parser.add_argument('--ymin',
                        type=float,
                        help='Min of plot.',
                        default=None)
    parser.add_argument('--ymax',
                        type=float,
                        help='Max of plot.',
                        default=None)
    args = parser.parse_args()

    if args.invert:
        op = lambda x: -x
    else:
        op = None

    print("Bhlight: Plotting {} ".format(args.qname),
          "in simulation directory {}.".format(args.dumpdir))
    plot_q_of_t(
        qname = args.qname,
        path = args.dumpdir,
        op = op,
        save = True,
        display = False,
        ymin = args.ymin,
        ymax = args.ymax,
        log = args.log
    )
    print("Done!")    
