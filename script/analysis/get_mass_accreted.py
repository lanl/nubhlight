#!/usr/bin/env python
# Author: Jonah Miller (jonahm@lanl.gov)

import numpy as np
from scipy import integrate
from units import cgs
import hdf5_to_dict as io

def get_m_accreted(path):
    "Given dump dir, get mass accreted, in solar masses"
    dumps = io.get_dumps_reduced(path,True)
    if len(dumps) == 0:
        dumps = io.get_dumps_reduced(path,False)
    if len(dumps) == 0:
        raise ValueError("No dumps to load!")
    hdr = io.load_hdr(dumps[0])
    diag = io.load_diag(path,hdr,timers=False)

    m_accreted = integrate.trapz(diag['mdot'],x=diag['t'])
    return np.abs(m_accreted*hdr['M_unit']/cgs['MSOLAR'])

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Get total mass accreted")
    parser.add_argument("dumpdirs",type=str,nargs="+",
                        help="directories to analyze")
    args = parser.parse_args()
    for d in args.dumpdirs:
        print(d,":",get_m_accreted(d))
