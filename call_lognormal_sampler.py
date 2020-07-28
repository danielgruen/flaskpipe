
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
parser.add_argument('--lognormal_mask_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/mask_sdss_sum.fits.gz")
parser.add_argument('--bin_edges_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/bin_edges.npy")
parser.add_argument('--galaxy_catalog_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits")
parser.add_argument('--galaxy_mask_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits")
parser.add_argument('--gamma_ray_map_unmasked_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/flux_9years_foresub")
parser.add_argument('--gamma_ray_mask_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/masks")
parser.add_argument('--exposure_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/expos_9years_binned")
options = parser.parse_args()
nside = int(options.nside)
nside_for_galaxy_mask_map = int(options.nside_for_galaxy_mask_map)
lognormal_mask_filename = options.lognormal_mask_filename
bin_edges_filename = options.bin_edges_filename
galaxy_catalog_filename = options.galaxy_catalog_filename
galaxy_mask_filename = options.galaxy_mask_filename
gamma_ray_map_unmasked_directory = options.gamma_ray_map_unmasked_directory
gamma_ray_mask_map_directory = options.gamma_ray_mask_map_directory
exposure_map_directory = options.exposure_map_directory
print("nside: {0}".format(nside))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("lognormal_mask_filename: {0}".format(lognormal_mask_filename))
print("bin_edges_filename: {0}".format(bin_edges_filename))
print("galaxy_catalog_filename: {0}".format(galaxy_catalog_filename))
print("galaxy_mask_filename: {0}".format(galaxy_mask_filename))
print("gamma_ray_map_unmasked_directory: {0}".format(gamma_ray_map_unmasked_directory))
print("gamma_ray_mask_map_directory: {0}".format(gamma_ray_mask_map_directory))
print("exposure_map_directory: {0}".format(exposure_map_directory))

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
    os.system("bsub -W 50 -R rhel60 -o {0}/lognormal_sampler.txt python {0}/lognormal_sampler.py --nside {1} --nside_for_galaxy_mask_map {2} --lognormal_mask_filename {3} --bin_edges_filename {4} --galaxy_catalog_filename {5} --galaxy_mask_filename {6} --gamma_ray_map_unmasked_directory {7} --gamma_ray_mask_map_directory {8}".format(directory, nside, lognormal_mask_filename, bin_edges_filename, galaxy_catalog_filename, galaxy_mask_filename, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory))


