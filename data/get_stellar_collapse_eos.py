#!/usr/bin/env python

"""get_eos_files.py
Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
Time-stamp: <2018-02-27 09:46:23 (jonahm)>

This is a quick and dirty script to download EOS files
from the Stellar Collapse Website.

Uses wget rather than Python's built in urllib
"""

from __future__ import print_function
from subprocess import call
from os import path
import sys, re, glob

URLS = {
    "https://stellarcollapse.org/EOS/"
    : ["LS180_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2"
       , "LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2"
       , "LS375_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2"
       , "KDE0v1_3335_rho391_temp163_ye66.h5.bz2"
       , "KDE0v1_0000_rho391_temp163_ye66.h5.bz2"
       , "LNS_3335_rho391_temp163_ye66.h5.bz2"
       , "LNS_0000_rho391_temp163_ye66.h5.bz2"
       , "LS220_3335_rho391_temp163_ye66.h5.bz2"
       , "LS220_0000_rho391_temp163_ye66.h5.bz2"
       , "LS220star_3335_rho391_temp163_ye66.h5.bz2"
       , "LS220star_0000_rho391_temp163_ye66.h5.bz2"
       , "NRAPR_3335_rho391_temp163_ye66.h5.bz2"
       , "NRAPR_0000_rho391_temp163_ye66.h5.bz2"
       , "SKRA_3335_rho391_temp163_ye66.h5.bz2"
       , "SKRA_0000_rho391_temp163_ye66.h5.bz2"
       , "SkT1_3335_rho391_temp163_ye66.h5.bz2"
       , "SkT1_0000_rho391_temp163_ye66.h5.bz2"
       , "SkT1star_3335_rho391_temp163_ye66.h5.bz2"
       , "SkT1star_0000_rho391_temp163_ye66.h5.bz2"
       , "Skxs20_3335_rho391_temp163_ye66.h5.bz2"
       , "Skxs20_0000_rho391_temp163_ye66.h5.bz2"
       , "SLy4_3335_rho391_temp163_ye66.h5.bz2"
       , "SLy4_0000_rho391_temp163_ye66.h5.bz2"
       , "SQMC700_3335_rho391_temp163_ye66.h5.bz2"
       , "SQMC700_0000_rho391_temp163_ye66.h5.bz2"
    ],
    "https://stellarcollapse.org/~evanoc/"
    : ["LS220_240r_140t_50y_analmu_20120628_SVNr28.h5.bz2"
       , "HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5.bz2"
       , "HShen_HyperonEOS_rho220_temp180_ye65_version_1.1_20131007.h5.bz2"
       , "GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5.bz2"
       , "GShenFSU_1.7EOS_rho280_temp180_ye52_version_1.1_20120817.h5.bz2"
       , "GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5.bz2"
       , "Hempel_TMAEOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "Hempel_TM1EOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "Hempel_FSGEOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "Hempel_NL3EOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "Hempel_IUFEOS_rho234_temp180_ye60_version_1.1_20140129.h5.bz2"
       , "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"
       , "BHB_lEOS_rho234_temp180_ye60_version_1.02_20140422.h5.bz2"
       , "BHB_lpEOS_rho234_temp180_ye60_version_1.02_20140422.h5.bz2"
    ]
}
       
def get_url(base,name):
    return "{}".format(path.join(base,name)).lstrip().rstrip()

def get_current_files():
    "Check what files are currently downloaded."
    h5files = [path.basename(f).rstrip('.h5') \
               for f in glob.glob("*.h5")]
    bz2files = [path.basename(f).rstrip('.h5.bz2') \
                for f  in glob.glob("*.h5.bz2")]
    all_files = h5files + bz2files
    return h5files,bz2files, all_files

def file_is_downloaded(name):
    h5files,bz2files,all_files = get_current_files()
    my_name = name.rstrip('.h5').rstrip('.h5.bz2')
    return my_name in all_files

def file_is_decompressed(name):
    h5files,bz2files,all_files = get_current_files()
    my_name = name.rstrip('.h5').rstrip('.h5.bz2')
    return my_name in h5files

def download(base,name):
    "Download a single EOS and decompress it."
    print("-------------------------------")
    print("Downloading file {}".format(name))
    url = get_url(base,name)
    if file_is_downloaded(name):
        print("\tFile already downloaded!")
    else:
        print("\tDownloading via wget.")
        call(["wget",url])
        print("\tDownload complete!")
    if file_is_decompressed(name):
        print("\tFile already decompressed.")
    else:
        print("\tAttempting to decompress")
        try:
            call(['bzip2','-d',name])
        except:
            print("\tDecompression failed.")
    print("\tDone.")
    print("-------------------------------")

def download_pattern(pattern):
    "Match to a pattern and download that EOS."
    regexp = re.compile(pattern)
    for base, namelist in URLS.items():
        for name in namelist:
            if regexp.search(name):
                download(base,name)

def download_all():
    "Download all URLS."
    for base,namelist in URLS.items():
        for name in namelist:
            download(base,name)
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description = ("Download EOS tables to be "
                       +"read by stellar collapse reader."),
        epilog = "Default is do nothing.")
    parser.add_argument('-a','--all',
                       dest='all',
                       action='store_true',
                       help='Download all files. Takes precedence over -p.')
    parser.add_argument('-p',
                        dest='pattern',
                        type=str,
                        default=None,
                        help='Download EOSs whose name match the pattern.')
    args = parser.parse_args()

    if args.all:
        print("Downloading all files.")
        download_all()
        print("Done!")
    elif args.pattern is not None:
        print("Downloading files that match pattern {}".format(args.pattern))
        download_pattern(args.pattern)
        print("Done!")
    else:
        parser.print_help(sys.stderr)
