
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

np.load(bin_edges_filename)
zmins = bin_edges[0:5]
zmaxs = bin_edges[1:6]
Emins = bin_edges[6,15]
Emaxs = bin_edges[7,16]

base_nside = 256
base_number_of_pixels = 786432.0
ratio_of_nside_to_base_nside = nside/base_nside
number_of_powers_above_base = np.log(ratio_of_nside_to_base_nside)/np.log(2.0)
number_of_pixels = base_number_of_pixels * (4.0 ** number_of_powers_above_base)
number_of_pixels_int = int(number_of_pixels)

np.random.seed(12345)





def pick_random_galon_and_galat_inside_pixel_of_known_index_and_location(galon_of_pixel, galat_of_pixel, theta_of_pixel, pix_number, pixel_half_sqrt_angular_area):
    pixel_half_sqrt_angular_area = np.sqrt((4 * np.pi)/number_of_pixels)/2.0
    expansive_angular_increment = pixel_half_sqrt_angular_area * 3.0
    galon_min = galon_of_pixel - expansive_angular_increment/np.sin(theta_of_pixel)
    galon_max = galon_of_pixel + expansive_angular_increment/np.sin(theta_of_pixel)
    galat_min = galat_of_pixel - expansive_angular_increment
    galat_max = galat_of_pixel + expansive_angular_increment
    while True:
        galon_of_pixel_perturbed = np.random.uniform(galon_min, galon_max)
        singalat_of_pixel_perturbed = np.random.uniform(np.sin(galat_min), np.sin(galat_max))
        galat_of_pixel_perturbed = np.arcsin(singalat_of_pixel_perturbed)

        theta_of_pixel_perturbed = -galat_of_pixel_perturbed + np.pi/2.0
        phi_of_pixel_perturbed = galon_of_pixel_perturbed
        potential_pix_number = hp.pixelfunc.ang2pix(nside=nside, theta = theta_of_pixel_perturbed, phi = phi_of_pixel_perturbed)
        if potential_pix_number == pix_number:
            return galon_of_pixel_perturbed, galat_of_pixel_perturbed






def get_unmasked_pixels_generic_given_generic_mask_map(generic_mask_map):
    unmasked_pixels_generic = []
    for index in list(range(0,number_of_pixels_int)):
        if generic_mask_map[index] > 0.0:
            unmasked_pixels_generic.append(index)
    print(len(unmasked_pixels_generic))
    return unmasked_pixels_generic





def get_nonzero_pixels_generic_given_generic_map(generic_map):
    nonzero_pixels_generic = []
    for index in list(range(0,number_of_pixels_int)):
        if generic_map[index] != 0.0:
            nonzero_pixels_generic.append(index)
    print(len(nonzero_pixels_generic))
    return nonzero_pixels_generic






def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension







def get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_or_unmasked_pixels_generic, generic_map_or_generic_mask_map, for_galaxies=False):
    pixel_half_sqrt_angular_area = np.sqrt((4 * np.pi)/number_of_pixels)/2.0
    thetas_of_pixels, phis_of_pixels = hp.pixelfunc.pix2ang(nside,nonzero_or_unmasked_pixels_generic)

    galons_of_pixels = phis_of_pixels
    galats_of_pixels = -thetas_of_pixels + np.pi/2.0

    galons_of_pixels_perturbed = []
    galats_of_pixels_perturbed = []
    weights_of_pixels_perturbed = []

    for index, galon_of_pixel, galat_of_pixel, theta_of_pixel, pix_number in zip(list(range(0, len(nonzero_or_unmasked_pixels_generic))), galons_of_pixels, galats_of_pixels, thetas_of_pixels, nonzero_or_unmasked_pixels_generic):
            if index % 10000 == 0:
                pass
            try:

                galon_of_pixel_perturbed, galat_of_pixel_perturbed = pick_random_galon_and_galat_inside_pixel_of_known_index_and_location(galon_of_pixel, galat_of_pixel, theta_of_pixel, pix_number, pixel_half_sqrt_angular_area)
                weight_of_pixel_perturbed = generic_map_or_generic_mask_map[pix_number]
                if weight_of_pixel_perturbed > 1 and for_galaxies:
                    for w in list(range(0,weight_of_pixel_perturbed)):
                        galon_of_pixel_perturbed, galat_of_pixel_perturbed = pick_random_galon_and_galat_inside_pixel_of_known_index_and_location(galon_of_pixel, galat_of_pixel, theta_of_pixel, pix_number, pixel_half_sqrt_angular_area)
                        galons_of_pixels_perturbed.append(galon_of_pixel_perturbed)
                        galats_of_pixels_perturbed.append(galat_of_pixel_perturbed)
                        weights_of_pixels_perturbed.append(weight_of_pixel_perturbed)                        
                else:
                    galons_of_pixels_perturbed.append(galon_of_pixel_perturbed)
                    galats_of_pixels_perturbed.append(galat_of_pixel_perturbed)
                    weights_of_pixels_perturbed.append(weight_of_pixel_perturbed)
            except:
                pass

    galons_of_pixels_perturbed = np.array(galons_of_pixels_perturbed)
    galats_of_pixels_perturbed = np.array(galats_of_pixels_perturbed)
    weights_of_pixels_perturbed = np.array(weights_of_pixels_perturbed)

    phis_gal = copy.deepcopy(galons_of_pixels_perturbed)
    thetas_gal = -copy.deepcopy(galats_of_pixels_perturbed) + np.pi/2.0
    thetas_equat, phis_equat = hp.rotator.Rotator(coord='gc')(thetas_gal, phis_gal)
    ras_of_pixels_perturbed = copy.deepcopy(phis_equat)
    decs_of_pixels_perturbed = -copy.deepcopy(thetas_equat) + np.pi/2.0

    return ras_of_pixels_perturbed, decs_of_pixels_perturbed, weights_of_pixels_perturbed






def generate_lognormal_galaxy_catalogues(zmin, zmax, hdu, number_of_unmasked_pixels, lognormal_mask, core_directory):
    ras = hdu[1].data['RA']
    zs = hdu[1].data['ZREDMAGIC']
    delete_list = []
    for i, z in enumerate(zs):
        if z < zmin or z > zmax:
            delete_list.append(i)

    ras = np.delete(ras, delete_list)
    number_of_galaxies = len(ras)

    number_of_galaxies_per_unmasked_pixel = number_of_galaxies/number_of_unmasked_pixels    
    

    filename = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/output_round0/sdss_12345/map-{0}.fits".format(fz_delineation)
    smooth_galaxy_overdensity_map = hp.read_map(filename)
    unsampled_galaxy_map = number_of_galaxies_per_unmasked_pixel * (smooth_galaxy_overdensity_map + 1.0)

    sampled_galaxy_map = np.random.poisson(unsampled_galaxy_map).astype(np.float)
    sampled_galaxy_map = sampled_galaxy_map * lognormal_mask
    nonzero_pixels_galaxy = get_nonzero_pixels_generic_given_generic_map(sampled_galaxy_map)
    ras, decs, weights = get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_pixels_galaxy, sampled_galaxy_map, for_galaxies)

    np.save("{0}/sampled_ras_{1}.npy".format(core_directory, fz_delineation), ras)
    np.save("{0}/sampled_decs_{1}.npy".format(core_directory, fz_delineation), decs)







def generate_lognormal_gamma_ray_catalogues(Emin, Emax, fz_delineation, lognormal_mask, core_directory):
    
    gamma_ray_map_unmasked = hp.read_map("/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/flux_9years_foresub/flux_9years_{0}_{1}_MeV_C.fits".format(Emin,Emax))
    gamma_ray_mask_map_high_res = hp.read_map("/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/masks/mask_GP30.0_sources_variable_FL8Y_incl_3FHL_incl_{0}_{1}_MeV_hpx_ord10.fits".format(Emin,Emax))
    gamma_ray_map_high_res = gamma_ray_map_unmasked * gamma_ray_mask_map_high_res



    gamma_ray_map = hp.ud_grade(gamma_ray_map_high_res, nside)
    exposure_map = hp.read_map("/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/expos_9years_binned/expos_9years_{0}_{1}_MeV_macro.fits".format(Emin,Emax))    
    photon_map = gamma_ray_map * exposure_map
    gamma_ray_mask_map = hp.ud_grade(gamma_ray_mask_map_high_res, nside)



    unmasked_pixels_gamma_ray = get_unmasked_pixels_generic_given_generic_mask_map(gamma_ray_mask_map)
    number_of_unmasked_pixels = len(unmasked_pixels_gamma_ray)
    number_of_photons_per_unmasked_pixel = np.sum(photon_map)/number_of_unmasked_pixels
    
    

    filename = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/output_round0/sdss_12345/map-{0}.fits".format(fz_delineation)
    smooth_photon_overdensity_map = hp.read_map(filename)
    unsampled_photon_map = number_of_photons_per_unmasked_pixel * (smooth_photon_overdensity_map + 1.0)

    sampled_photon_map = np.random.poisson(unsampled_photon_map).astype(np.float)
    sampled_gamma_ray_map = sampled_photon_map / exposure_map
    
    hp.write_map("{0}/sampled_gamma_ray_map".format(core_directory, fz_delineation), sampled_gamma_ray_map)








core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()
lognormal_mask = hp.read_map(lognormal_mask)

for Emin, Emax in zip(Emins, Emaxs):
    os.system("mkdir {0}/new_foreground_subtraction_{1}_{2}_MeV".format(core_directory, Emin, Emax))
    
run_directories = []
for Emin, Emax in zip(Emins, Emaxs):
    directory = "{0}/new_foreground_subtraction_{1}_{2}_MeV".format(core_directory, Emin, Emax)
    for zmin, zmax in zip(zmins, zmaxs):
        inner_directory = "mkdir {0}/z_{1}_{2}".format(directory, zmin, zmax)
        os.system("mkdir {}".format(inner_directory))
    os.system("cp {0}/call_cross_correlations.py {0}/cross_correlations.py {1}".format(core_directory,directory))
    run_directories.append(directory)




RM_galaxy_mask_map = np.zeros(201326592)
filename = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits"
hdu = fits.open(filename)
for index in list(range(0,len(hdu[1].data["HPIX"]))):
    unmasked_pixel = hdu[1].data["HPIX"][index]
    RM_galaxy_mask_map[unmasked_pixel] = hdu[1].data["FRACGOOD"][index]

unmasked_pixels = []
for index in list(range(0,201326592)):
    if RM_galaxy_mask_map_low_res[index] > 0.0:
        unmasked_pixels.append(index)
unmasked_pixels = np.array(unmasked_pixels)
number_of_unmasked_pixels = len(unmasked_pixels)
print(number_of_unmasked_pixels)







filename = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits"
hdu = fits.open(filename)

#zmins = [0.15, 0.35, 0.5, 0.65, 0.8]
#zmaxs = [0.35, 0.5, 0.65, 0.8, 0.9]
fz_delineations = ["f1z1", "f2z2", "f3z3", "f4z4", "f5z5"]
for zmin, zmax, fz_delineation in zip(zmins, zmaxs, fz_delineations):
    generate_lognormal_galaxy_catalogues(zmin, zmax, hdu, number_of_unmasked_pixels)

#Emins = [631.0, 1202.3, 2290.9, 4786.3, 9120.1, 17378.0, 36307.8, 69183.1, 131825.7]
#Emaxs = [1202.3, 2290.9, 4786.3, 9120.1, 17378.0, 36307.8, 69183.1, 131825.7, 1000000.0]
fz_delineations = ["f6z1", "f6z2", "f6z3", "f6z4", "f6z5", "f6z6", "f6z7", "f6z8", "f6z9"]
for Emin, Emax, fz_delineation in zip(Emins, Emaxs, fz_delineations):
    generate_lognormal_gamma_ray_catalogues(Emin, Emax, fz_delineation)



for run_directory in run_directories:
    os.system("bsub -o {0}/call_cross_correlations.txt python {0}/call_cross_correlations.py --nside {1} --zbin_edge1 {2} --zbin_edge2 {3} --zbin_edge3 {4} --zbin_edge4 {5} --zbin_edge5 {6} --zbin_edge6 {7} --data_type {8}".format(run_directory,nside,bin_edges[0],bin_edges[1],bin_edges[2],bin_edges[3],bin_edges[4],bin_edges[5],"log_realiz"))
