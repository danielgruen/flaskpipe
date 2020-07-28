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
parser.add_argument('--nside_for_galaxy_mask_map')
parser.add_argument('--bin_edges_filename')
parser.add_argument('--data_type')
parser.add_argument('--Emin')
parser.add_argument('--Emax')
parser.add_argument('--galaxy_catalog_filename')
parser.add_argument('--galaxy_mask_filename')
parser.add_argument('--gamma_ray_map_unmasked_directory')
parser.add_argument('--gamma_ray_mask_map_directory')
options = parser.parse_args()
nside = int(options.nside)
nside_for_galaxy_mask_map = int(options.nside_for_galaxy_mask_map)
bin_edges_filename = options.bin_edges_filename
data_type = options.data_type
Emin = options.Emin
Emax = options.Emax
galaxy_catalog_filename = options.galaxy_catalog_filename
galaxy_mask_filename = options.galaxy_mask_filename
gamma_ray_map_unmasked_directory = options.gamma_ray_map_unmasked_directory
gamma_ray_mask_map_directory = options.gamma_ray_mask_map_directory
print("nside: {0}".format(nside))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("bin_edges_filename: {0}".format(bin_edges_filename))
print("data_type: {0}".format(data_type))
print("Emin: {0}".format(Emin))
print("Emax: {0}".format(Emax))
print("galaxy_catalog_filename: {0}".format(galaxy_catalog_filename))
print("galaxy_mask_filename: {0}".format(galaxy_mask_filename))
print("gamma_ray_map_unmasked_directory: {0}".format(gamma_ray_map_unmasked_directory))
print("gamma_ray_mask_map_directory: {0}".format(gamma_ray_mask_map_directory))




def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension

core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

print("bsub -W 50 -R rhel60 -o {0}/cross_correlations.txt python {0}/cross_correlations.py --nside {1} --nside_for_galaxy_mask_map {2} --bin_edges_filename {3} --data_type {4} --Emin {5} --Emax {6} --galaxy_catalog_filename {7} --galaxy_mask_filename {8} --gamma_ray_map_unmasked_directory {9} --gamma_ray_mask_map_directory {10}".format(core_directory, nside, nside_for_galaxy_mask_map, bin_edges_filename, data_type, Emin, Emax, galaxy_catalog_filename, galaxy_mask_filename, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory))

os.system("bsub -W 50 -R rhel60 -o {0}/cross_correlations.txt python {0}/cross_correlations.py --nside {1} --nside_for_galaxy_mask_map {2} --bin_edges_filename {3} --data_type {4} --Emin {5} --Emax {6} --galaxy_catalog_filename {7} --galaxy_mask_filename {8} --gamma_ray_map_unmasked_directory {9} --gamma_ray_mask_map_directory {10}".format(core_directory, nside, nside_for_galaxy_mask_map, bin_edges_filename, data_type, Emin, Emax, galaxy_catalog_filename, galaxy_mask_filename, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory))


