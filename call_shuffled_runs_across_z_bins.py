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
options = parser.parse_args()
nside = int(options.nside)
bin_edges_filename = options.bin_edges_filename
print("nside: {0}".format(nside))
print("bin_edges_filename: {0}".format(bin_edges_filename))

bin_edges = np.load(bin_edges_filename)
zmins = bin_edges[0:5]
zmaxs = bin_edges[1:6]

def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension

core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

for zmin, zmax in zip(zmins, zmaxs):
    folder_name = "z_{0}_{1}/shuffled_new_foreground_subtraction".format(zmin,zmax)
    os.system("mkdir {0}".format(folder_name))
    os.system("cp call_shuffled_runs.py {0}".format(folder_name))
    os.system("cp shuffled_run.py {0}".format(folder_name))
    os.system("bsub -W 10 -R rhel60 -o {0}/{1}/call_shuffled_runs.txt python {0}/{1}/call_shuffled_runs.py --nside {2} --zmin {3} --zmax {4}".format(core_directory, folder_name, nside, zmin, zmax))
    print("bsub -W 10 -R rhel60 -o {0}/{1}/call_shuffled_runs.txt python {0}/{1}/call_shuffled_runs.py --nside {2} --zmin {3} --zmax {4}".format(core_directory, folder_name, nside, zmin, zmax))

