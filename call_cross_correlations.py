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
parser.add_argument('--bin_edges_filename')
parser.add_argument('--data_type')
options = parser.parse_args()
nside = int(options.nside)
bin_edges_filename = options.bin_edges_filename
data_type = float(options.data_type)
print("nside: {0}".format(nside))
print("bin_edges_filename: {0}".format(bin_edges_filename))
print("data_type: {0}".format(data_type))

np.load(bin_edges_filename)
zmins = bin_edges[0:5]
zmaxs = bin_edges[1:6]
Emins = bin_edges[6,15]
Emaxs = bin_edges[7,16]




def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension

core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

os.system("bsub -W 79 -R rhel60 -o {0}/cross_correlations.txt python {0}/cross_correlations.py --nside {1} --bin_edges_filename {2} --data_type {3}".format(core_directory, nside, bin_edges_filename, data_type))


