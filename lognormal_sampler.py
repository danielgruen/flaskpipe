
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
parser.add_argument('--lognormal_mask_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/mask_sdss_sum.fits.gz")
parser.add_argument('--bin_edges_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/bin_edges.npy")
parser.add_argument('--galaxy_catalog_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits")
parser.add_argument('--galaxy_mask_filename', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits")
parser.add_argument('--gamma_ray_map_unmasked_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/flux_9years_foresub")
parser.add_argument('--gamma_ray_mask_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/masks")
parser.add_argument('--exposure_map_directory', default = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/expos_9years_binned")
options = parser.parse_args()
nside = int(options.nside)
lognormal_mask_filename = options.lognormal_mask_filename
bin_edges_filename = options.bin_edges_filename
galaxy_catalog_filename = options.galaxy_catalog_filename
galaxy_mask_filename = options.galaxy_mask_filename
gamma_ray_map_unmasked_directory = options.gamma_ray_map_unmasked_directory
gamma_ray_mask_map_directory = options.gamma_ray_mask_map_directory
exposure_map_directory = options.exposure_map_directory
print("nside: {0}".format(nside))
print("lognormal_mask_filename: {0}".format(lognormal_mask_filename))
print("bin_edges_filename: {0}".format(bin_edges_filename))
print("galaxy_catalog_filename: {0}".format(galaxy_catalog_filename))
print("galaxy_mask_filename: {0}".format(galaxy_mask_filename))
print("gamma_ray_map_unmasked_directory: {0}".format(gamma_ray_map_unmasked_directory))
print("gamma_ray_mask_map_directory: {0}".format(gamma_ray_mask_map_directory))
print("exposure_map_directory: {0}".format(exposure_map_directory))

np.load(bin_edges_filename)
zmins = bin_edges[0:5]
zmaxs = bin_edges[1:6]
Emins = bin_edges[6,15]
Emaxs = bin_edges[7,16]

number_of_pixels_int = hp.nside2npix(nside)
number_of_pixels = float(number_of_pixels_int)

np.random.seed(12345)





def get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_or_unmasked_pixels_generic, generic_map_or_generic_mask_map, for_sampled_galaxies=False):

    if for_sampled_galaxies:
        number_of_random_points = np.sum(generic_map_or_generic_mask_map)
        pixels_for_all_random_points = []
        counter = 0
        while True:
            pixels_with_more_than_counter_galaxies = np.where(generic_map_or_generic_mask_map > counter)[0].tolist()
            if len(pixels_with_more_than_counter_galaxies) == 0:
                break
            pixels_for_all_random_points = pixels_for_all_random_points + pixels_with_more_than_counter_galaxies
            counter = counter + 1 
        pixels_for_all_random_points = np.array(pixels_for_all_random_points)
    else:
        number_of_random_points = len(nonzero_or_unmasked_pixels_generic)
        pixels_for_all_random_points = copy.deepcopy(nonzero_or_unmasked_pixels_generic)



    pixel_half_sqrt_angular_area = np.sqrt((4 * np.pi)/number_of_pixels)/2.0
    expansive_angular_increment = pixel_half_sqrt_angular_area * 3.0
    thetas_of_pixels, phis_of_pixels = hp.pixelfunc.pix2ang(nside, pixels_for_all_random_points)
    galons_of_pixels = phis_of_pixels
    galats_of_pixels = -thetas_of_pixels + np.pi/2.0
    galons_of_pixels_perturbed = []
    galats_of_pixels_perturbed = []
    weights_of_pixels_perturbed = []


    pixels_with_locations_not_found_all_random_points = np.ones(number_of_random_points, dtype=bool)
    while True:        
        thetas_with_locations_not_found = thetas_of_pixels[pixels_with_locations_not_found_all_random_points]
        galons_with_locations_not_found = galons_of_pixels[pixels_with_locations_not_found_all_random_points]
        galats_with_locations_not_found = galats_of_pixels[pixels_with_locations_not_found_all_random_points]

        galon_mins = galons_with_locations_not_found - expansive_angular_increment/np.sin(thetas_with_locations_not_found)
        galon_maxs = galons_with_locations_not_found + expansive_angular_increment/np.sin(thetas_with_locations_not_found)
        galat_mins = galats_with_locations_not_found - expansive_angular_increment
        galat_maxs = galats_with_locations_not_found + expansive_angular_increment

        galons_of_pixels_perturbed_with_locations_not_found = np.random.uniform(galon_mins, galon_maxs)
        singalats_of_pixels_perturbed = np.random.uniform(np.sin(galat_mins), np.sin(galat_maxs))
        galats_of_pixels_perturbed_with_locations_not_found = np.arcsin(singalats_of_pixels_perturbed)

        thetas_of_pixels_perturbed = -galats_of_pixels_perturbed_with_locations_not_found + np.pi/2.0
        phis_of_pixels_perturbed = copy.deepcopy(galons_of_pixels_perturbed_with_locations_not_found)
        potential_pix_numbers = hp.pixelfunc.ang2pix(nside=nside, theta = thetas_of_pixels_perturbed, phi = phis_of_pixels_perturbed)
        successful_match_indices = np.where(potential_pix_numbers == pixels_with_locations_not_found_all_random_points)[0]
        failed_match_indices = np.where(potential_pix_numbers != pixels_with_locations_not_found_all_random_points)[0]

        galons_of_pixels_perturbed = galons_of_pixels_perturbed + galons_of_pixels_perturbed_with_locations_not_found[successful_match_indices].tolist()
        galats_of_pixels_perturbed = galats_of_pixels_perturbed + galats_of_pixels_perturbed_with_locations_not_found[successful_match_indices].tolist()
        pixels_with_locations_just_found_all_random_points = pixels_with_locations_not_found_all_random_points[successful_match_indices]
        weights_of_pixels_perturbed = weights_of_pixels_perturbed + generic_map_or_generic_mask_map[pixels_with_locations_just_found_all_random_points].tolist()
        pixels_with_locations_not_found_all_random_points = pixels_with_locations_not_found_all_random_points[failed_match_indices]

    galons_of_pixels_perturbed = np.array(galons_of_pixels_perturbed)
    galats_of_pixels_perturbed = np.array(galats_of_pixels_perturbed)
    weights_of_pixels_perturbed = np.array(weights_of_pixels_perturbed)
    if for_sampled_galaxies:
        weights_of_pixels_perturbed = np.ones(number_of_random_points)

    phis_gal = copy.deepcopy(galons_of_pixels_perturbed)
    thetas_gal = -copy.deepcopy(galats_of_pixels_perturbed) + np.pi/2.0
    thetas_equat, phis_equat = hp.rotator.Rotator(coord='gc')(thetas_gal, phis_gal)
    ras_of_pixels_perturbed = copy.deepcopy(phis_equat)
    decs_of_pixels_perturbed = -copy.deepcopy(thetas_equat) + np.pi/2.0

    return ras_of_pixels_perturbed, decs_of_pixels_perturbed, weights_of_pixels_perturbed



def generate_lognormal_galaxy_catalogues(zmin, zmax, fz_delineation, hdu_data_table, sum_of_pixel_mask_values, lognormal_mask, core_directory):

    ras = hdu_data_table['RA']
    zs = hdu_data_table['ZREDMAGIC']
    new_ras_indices = np.where((zs > zmin) & (zs < zmax))[0]
    number_of_galaxies = len(new_ras_indices)
    number_of_galaxies_per_unmasked_pixel = number_of_galaxies/sum_of_pixel_mask_values    
    

    smooth_galaxy_overdensity_map = hp.read_map("{0}/map-{1}.fits".format(core_directory,fz_delineation))
    smooth_galaxy_overdensity_map = hp.ud_grade(smooth_galaxy_overdensity_map, nside)
    if np.any(smooth_galaxy_overdensity_map < -1.0):
        warnings.warn("Warning: smooth_galaxy_overdensity_map has values < -1.0. Setting these to -1.0 to prevent issues with the Poisson sampling")
        smooth_galaxy_overdensity_map = np.where(smooth_galaxy_overdensity_map >= -1.0, smooth_galaxy_overdensity_map, -1.0)
    unsampled_galaxy_map = number_of_galaxies_per_unmasked_pixel * (smooth_galaxy_overdensity_map + 1.0) * lognormal_mask

    sampled_galaxy_map = np.random.poisson(unsampled_galaxy_map).astype(np.float)
    sampled_galaxy_map = sampled_galaxy_map
    nonzero_pixels_galaxy = np.where(sampled_galaxy_map != 0.0)[0]
    ras, decs, weights = get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_pixels_galaxy, sampled_galaxy_map, for_sampled_galaxies=True)

    np.save("{0}/sampled_ras_z_{1}_{2}.npy".format(core_directory, zmin, zmax), ras)
    np.save("{0}/sampled_decs_z_{1}_{2}.npy".format(core_directory, zmin, zmax), decs)



def generate_lognormal_gamma_ray_catalogues(Emin, Emax, fz_delineation, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory, exposure_map_directory, lognormal_mask, core_directory, nside):
    
    gamma_ray_map_unmasked = hp.read_map("{0}/flux_9years_{1}_{2}_MeV_C.fits".format(gamma_ray_map_unmasked_directory, Emin,Emax)) 
    gamma_ray_mask_map_high_res = hp.read_map("{0}/mask_GP30.0_sources_variable_FL8Y_incl_3FHL_incl_{1}_{2}_MeV_hpx_ord10.fits".format(gamma_ray_mask_map_directory, Emin,Emax))
    gamma_ray_map_high_res = gamma_ray_map_unmasked * hp.ud_grade(gamma_ray_mask_map_high_res, hp.npix2nside(len(gamma_ray_map_unmasked)))
    gamma_ray_map = hp.ud_grade(gamma_ray_map_high_res, nside)
    exposure_map = hp.read_map("{0}/expos_9years_{1}_{2}_MeV_macro.fits".format(exposure_map_directory, Emin,Emax)) 
    exposure_map = hp.ud_grade(exposure_map, nside)   
    photon_map = gamma_ray_map * exposure_map
    gamma_ray_mask_map = hp.ud_grade(gamma_ray_mask_map_high_res, nside)
    sum_of_pixel_mask_values = np.sum(gamma_ray_mask_map)
    number_of_photons_per_unmasked_pixel = np.sum(photon_map)/sum_of_pixel_mask_values
    
    
    smooth_photon_overdensity_map = hp.read_map("{0}/map-{1}.fits".format(core_directory,fz_delineation))
    smooth_photon_overdensity_map = hp.ud_grade(smooth_photon_overdensity_map, nside)
    if np.any(smooth_photon_overdensity_map < -1.0):
        warnings.warn("Warning: smooth_photon_overdensity_map has values < -1.0. Setting these to -1.0 to prevent issues with the Poisson sampling")
    	smooth_photon_overdensity_map = np.where(smooth_photon_overdensity_map >= -1.0, smooth_photon_overdensity_map, -1.0)
    unsampled_photon_map = number_of_photons_per_unmasked_pixel * (smooth_photon_overdensity_map + 1.0) * lognormal_mask

    sampled_photon_map = np.random.poisson(unsampled_photon_map).astype(np.float)
    sampled_gamma_ray_map = sampled_photon_map / exposure_map
    
    hp.write_map("{0}/sampled_gamma_ray_map_{1}_{2}_MeV".format(core_directory, Emin, Emax), sampled_gamma_ray_map)



def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension





core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()
lognormal_mask = hp.read_map(lognormal_mask_filename)
lognormal_mask = hp.ud_grade(lognormal_mask, nside)

for Emin, Emax in zip(Emins, Emaxs):
    os.system("mkdir {0}/E_{1}_{2}_MeV".format(core_directory, Emin, Emax))
    
run_directories = []
for Emin, Emax in zip(Emins, Emaxs):
    directory = "{0}/E_{1}_{2}_MeV".format(core_directory, Emin, Emax)
    for zmin, zmax in zip(zmins, zmaxs):
        inner_directory = "mkdir {0}/z_{1}_{2}".format(directory, zmin, zmax)
        os.system("mkdir {}".format(inner_directory))
    os.system("cp {0}/call_cross_correlations.py {0}/cross_correlations.py {1}".format(core_directory,directory))
    run_directories.append(directory)




hdu = fits.open(galaxy_mask_filename)
RM_galaxy_mask_map = copy.deepcopy(hdu[1].data["FRACGOOD"])
RM_galaxy_mask_map = hp.ud_grade(RM_galaxy_mask_map, nside)
sum_of_pixel_mask_values = np.sum(RM_galaxy_mask_map)
#print(sum_of_pixel_mask_values)
hdu = fits.open(galaxy_catalog_filename)
hdu_data_table = hdu[1].data

#zmins = [0.15, 0.35, 0.5, 0.65, 0.8]
#zmaxs = [0.35, 0.5, 0.65, 0.8, 0.9]
fz_delineations = ["f1z1", "f2z2", "f3z3", "f4z4", "f5z5"]
for zmin, zmax, fz_delineation in zip(zmins, zmaxs, fz_delineations):
    generate_lognormal_galaxy_catalogues(zmin, zmax, fz_delineation, hdu_data_table, sum_of_pixel_mask_values, lognormal_mask, core_directory)




#Emins = [631.0, 1202.3, 2290.9, 4786.3, 9120.1, 17378.0, 36307.8, 69183.1, 131825.7]
#Emaxs = [1202.3, 2290.9, 4786.3, 9120.1, 17378.0, 36307.8, 69183.1, 131825.7, 1000000.0]
fz_delineations = ["f6z1", "f6z2", "f6z3", "f6z4", "f6z5", "f6z6", "f6z7", "f6z8", "f6z9"]
for Emin, Emax, fz_delineation in zip(Emins, Emaxs, fz_delineations):
    generate_lognormal_gamma_ray_catalogues(Emin, Emax, fz_delineation, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory, exposure_map_directory, lognormal_mask, core_directory, nside)




for run_directory, Emin, Emax in zip(run_directories, Emins, Emaxs):
    os.system("bsub -o {0}/call_cross_correlations.txt python {0}/call_cross_correlations.py --nside {1} --bin_edges_filename {2} --data_type {3} --Emin {4} --Emax {5} --galaxy_catalog_filename {6} --galaxy_mask_filename {7} --gamma_ray_map_unmasked_directory {8} --gamma_ray_mask_map_directory {9}".format(run_directory, nside, bin_edges_filename, "log_realiz", Emin, Emax, galaxy_catalog_filename, galaxy_mask_filename, gamma_ray_map_unmasked_directory, gamma_ray_mask_map_directory))
