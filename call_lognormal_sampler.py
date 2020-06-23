
from __future__ import print_function, division

import numpy as np
import pandas as pd
#from astropy.io import fits
import matplotlib.pyplot as plt
import os
import scipy.stats as st
import math
import glob

import lmfit
#import galsim
import piff

#from zernike import Zernike
import copy

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

import treecorr as tr
import healpy as hp
import fitsio
from astropy.io import fits
import fileinput
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--nside')
parser.add_argument('--lognormal_mask_filename')
parser.add_argument('--bin_edges_filename')
options = parser.parse_args()
nside = int(options.nside)
lognormal_mask_filename = options.lognormal_mask_filename
bin_edges_filename = options.bin_edges_filename
print("nside: {0}".format(nside))
print("lognormal_mask_filename: {0}".format(lognormal_mask_filename))
print("bin_edges_filename: {0}".format(bin_edges_filename))

def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension

core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

directories = glob.glob("{0}/output/*".format(core_directory))


for directory in directories:
    os.system("cp lognormal_sampler.py {0}".format(directory))
    os.system("cp call_cross_correlations.py {0}".format(directory))
    os.system("cp cross_correlations.py {0}".format(directory))
    os.system("bsub -W 79 -R rhel60 -o {0}/lognormal_sampler.txt python {0}/lognormal_sampler.py --nside {1} --lognormal_mask_filename {2} --bin_edges_filename {3}".format(directory, nside, lognormal_mask_filename, bin_edges_filename))


