from __future__ import print_function, division

# fix for DISPLAY variable issue
import matplotlib
matplotlib.use('Agg')

import numpy as np
from scipy import stats
from astropy.io import fits
import fitsio
import matplotlib.pyplot as plt
import os
import subprocess
import glob

import lmfit
import galsim
import piff

#from zernike import Zernike
import copy

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

import argparse
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser()
parser.add_argument('--nside')
parser.add_argument('--nside_for_galaxy_mask_map')
parser.add_argument('--bin_edges_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/bin_edges.npy")
parser.add_argument('--galaxy_catalog_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits")
parser.add_argument('--galaxy_mask_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits")
parser.add_argument('--gamma_ray_map_unmasked_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/flux_9years_foresub")
parser.add_argument('--gamma_ray_mask_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/masks")
options = parser.parse_args()
nside = int(options.nside)
nside_for_galaxy_mask_map = int(options.nside_for_galaxy_mask_map)
bin_edges_filename = options.bin_edges_filename
galaxy_catalog_filename = options.galaxy_catalog_filename
galaxy_mask_filename = options.galaxy_mask_filename
gamma_ray_map_unmasked_directory = options.gamma_ray_map_unmasked_directory
gamma_ray_mask_map_directory = options.gamma_ray_mask_map_directory
print("nside: {0}".format(nside))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("bin_edges_filename: {0}".format(bin_edges_filename))
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

def find_Emin_and_Emax_given_core_directory(core_directory):
    energy_range = core_directory.split("new_foreground_subtraction_")[-1]
    Emin = energy_range.split("_")[0]
    Emax = energy_range.split("_")[1]
    return Emin, Emax

core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

Emin, Emax = find_Emin_and_Emax_given_core_directory(core_directory)



for seed in list(range(12345,12395)):

        core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

        terminal_command = "bsub -W 79 -R rhel60 -o {0}/shuffled_run_seed_{1}.txt python {0}/shuffled_run.py --nside {2} --seed {1} --nside_for_galaxy_mask_map {3} --bin_edges_filename {4} --Emin {5} --Emax {6} --galaxy_catalog_filename {7} --galaxy_mask_filename {8} --gamma_ray_map_unmasked_directory {9} --gamma_ray_mask_map_directory {10}".format(core_directory, seed, nside, nside_for_galaxy_mask_map, bin_edges_filename, Emin, Emax, galaxy_catalog_filename, galaxy_mask_filename, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory)

        os.system(terminal_command)
        print(terminal_command)
