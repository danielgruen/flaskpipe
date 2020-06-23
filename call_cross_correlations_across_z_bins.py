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
parser.add_argument('--zbin_edge1')
parser.add_argument('--zbin_edge2')
parser.add_argument('--zbin_edge3')
parser.add_argument('--zbin_edge4')
parser.add_argument('--zbin_edge5')
parser.add_argument('--zbin_edge6')
parser.add_argument('--data_type')
options = parser.parse_args()
nside = int(options.nside)
zbin_edge1 = float(options.zbin_edge1)
zbin_edge2 = float(options.zbin_edge2)
zbin_edge3 = float(options.zbin_edge3)
zbin_edge4 = float(options.zbin_edge4)
zbin_edge5 = float(options.zbin_edge5)
zbin_edge6 = float(options.zbin_edge6)
data_type = float(options.data_type)
print("nside: {0}".format(nside))
print("zbin_edge1: {0}".format(zbin_edge1))
print("zbin_edge2: {0}".format(zbin_edge2))
print("zbin_edge3: {0}".format(zbin_edge3))
print("zbin_edge4: {0}".format(zbin_edge4))
print("zbin_edge5: {0}".format(zbin_edge5))
print("zbin_edge6: {0}".format(zbin_edge6))
print("data_type: {0}".format(data_type))

def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension

core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

zmins = [zbin_edge1, zbin_edge2, zbin_edge3, zbin_edge4, zbin_edge5]
zmaxs = [zbin_edge2, zbin_edge3, zbin_edge4, zbin_edge5, zbin_edge6]



os.system("bsub -W 79 -R rhel60 -o {0}/cross_correlations.txt python {0}/cross_correlations.py --nside {1} --zbin_edge1 {2} --zbin_edge2 {3} --zbin_edge3 {4} --zbin_edge4 {5} --zbin_edge5 {6} --zbin_edge6 {7} --data_type {8}".format(run_directory,nside,zbin_edge1,zbin_edge2,zbin_edge3,zbin_edge4,zbin_edge5,zbin_edge6,data_type))


