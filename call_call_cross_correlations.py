
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

#from zernike import Zernike
import copy
import sys

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

import treecorr as tr
import healpy as hp
import fitsio
from astropy.io import fits
import fileinput

import warnings
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--nside')
parser.add_argument('--nside_for_galaxy_mask_map')
parser.add_argument('--bin_edges_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/bin_edges.npy")
parser.add_argument('--galaxy_catalog_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits")
parser.add_argument('--galaxy_mask_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits")
parser.add_argument('--gamma_ray_map_unmasked_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/flux_9years_foresub")
parser.add_argument('--gamma_ray_mask_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/masks")
parser.add_argument('--exposure_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/expos_9years_binned")
options = parser.parse_args()
nside = int(options.nside)
nside_for_galaxy_mask_map = int(options.nside_for_galaxy_mask_map)
bin_edges_filename = options.bin_edges_filename
galaxy_catalog_filename = options.galaxy_catalog_filename
galaxy_mask_filename = options.galaxy_mask_filename
gamma_ray_map_unmasked_directory = options.gamma_ray_map_unmasked_directory
gamma_ray_mask_map_directory = options.gamma_ray_mask_map_directory
exposure_map_directory = options.exposure_map_directory
print("nside: {0}".format(nside))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("bin_edges_filename: {0}".format(bin_edges_filename))
print("galaxy_catalog_filename: {0}".format(galaxy_catalog_filename))
print("galaxy_mask_filename: {0}".format(galaxy_mask_filename))
print("gamma_ray_map_unmasked_directory: {0}".format(gamma_ray_map_unmasked_directory))
print("gamma_ray_mask_map_directory: {0}".format(gamma_ray_mask_map_directory))
print("exposure_map_directory: {0}".format(exposure_map_directory))

bin_edges = np.load(bin_edges_filename)
zmins = bin_edges[0:5]
zmaxs = bin_edges[1:6]
print(bin_edges)
Emins = bin_edges[6:15]
Emaxs = bin_edges[7:16]

number_of_pixels_int = hp.nside2npix(nside)
number_of_pixels = float(number_of_pixels_int)

np.random.seed(12345)





def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension





core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

for Emin, Emax in zip(Emins, Emaxs):
    print("mkdir {0}/E_{1}_{2}_MeV".format(core_directory, Emin, Emax))
    os.system("mkdir {0}/E_{1}_{2}_MeV".format(core_directory, Emin, Emax))
    
run_directories = []
for Emin, Emax in zip(Emins, Emaxs):
    directory = "{0}/E_{1}_{2}_MeV".format(core_directory, Emin, Emax)
    for zmin, zmax in zip(zmins, zmaxs):
        inner_directory = "{0}/z_{1}_{2}".format(directory, zmin, zmax)
        print("mkdir {}".format(inner_directory))
        os.system("mkdir {}".format(inner_directory))
    print("cp {0}/call_cross_correlations.py {0}/cross_correlations.py {1}".format(core_directory,directory))
    os.system("cp {0}/call_cross_correlations.py {0}/cross_correlations.py {1}".format(core_directory,directory))
    run_directories.append(directory)





for run_directory, Emin, Emax in zip(run_directories, Emins, Emaxs):
    os.system("bsub -o {0}/call_cross_correlations.txt python {0}/call_cross_correlations.py --nside {1} --bin_edges_filename {2} --data_type {3} --Emin {4} --Emax {5} --galaxy_catalog_filename {6} --galaxy_mask_filename {7} --gamma_ray_map_unmasked_directory {8} --gamma_ray_mask_map_directory {9}".format(run_directory, nside, bin_edges_filename, "log_realiz", Emin, Emax, fz_delineation, galaxy_catalog_filename, galaxy_mask_filename, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory))
