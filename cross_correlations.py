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
parser.add_argument('--Emin')
parser.add_argument('--Emax')
parser.add_argument('--gamma_catalog_filename')
parser.add_argument('--galaxy_mask_filename')
parser.add_argument('--gamma_ray_map_unmasked_directory')
parser.add_argument('--gamma_ray_mask_map_directory')
options = parser.parse_args()
nside = int(options.nside)
bin_edges_filename = float(options.bin_edges_filename)
data_type = options.data_type
Emin = options.Emin
Emax = options.Emax
gamma_catalog_filename = options.gamma_catalog_filename
galaxy_mask_filename = options.galaxy_mask_filename
gamma_ray_map_unmasked_directory = options.gamma_ray_map_unmasked_directory
gamma_ray_mask_map_directory = options.gamma_ray_mask_map_directory
print("nside: {0}".format(nside))
print("bin_edges_filename: {0}".format(bin_edges_filename))
print("data_type: {0}".format(data_type))
print("Emin: {0}".format(Emin))
print("Emax: {0}".format(Emax))
print("galaxy_catalog_filename: {0}".format(galaxy_catalog_filename))
print("galaxy_mask_filename: {0}".format(galaxy_mask_filename))
print("gamma_ray_map_unmasked_directory: {0}".format(gamma_ray_map_unmasked_directory))
print("gamma_ray_mask_map_directory: {0}".format(gamma_ray_mask_map_directory))

np.load(bin_edges_filename)
zmins = bin_edges[0:5]
zmaxs = bin_edges[1:6]

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



def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension





super_core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

    if data_type =="log_realiz":
        ras = np.load("{0}/../sampled_ras_z_{1}_{2}.npy".format(super_core_directory, zmin, zmax))
        decs = np.load("{0}/../sampled_decs_z_{1}_{2}.npy".format(super_core_directory, zmin, zmax))
        weights = np.ones(len(ras))
    else:
        hdu = fits.open("{0}/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits".format(galaxy_catalog_filename))

        ras = hdu[1].data['RA']
        decs = hdu[1].data['DEC']
        weights = hdu[1].data['weight']
        zs = hdu[1].data['ZREDMAGIC']
        new_ras_indices = np.where((zs > zmin) & (zs < zmax))[0]
                
        ras = ras[new_ras_indices]
        decs = decs[new_ras_indices]
        weights =  weights[new_ras_indices]

    cat_galaxy = tr.Catalog(ra=ras, dec=decs, w=weights, ra_units='degrees', dec_units='degrees')





for zmin, zmax in zip(zmins,zmaxs):
    core_directory = super_core_directory + "/z_{0}_{1}".format(zmin,zmax)
    if data_type =="log_realiz":
        hdu = fits.open("{0}/../../../mask_sdss_sum.fits.gz".format(super_core_directory))
    else:
        hdu = fits.open("{0}/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits".format(galaxy_mask_map_directory)
        hdu = fits.open(filename)
    RM_galaxy_mask_map = np.zeros(len(hdu[1].data["HPIX"]))
    for index in list(range(0,len(hdu[1].data["HPIX"]))):
        unmasked_pixel = hdu[1].data["HPIX"][index]
        RM_galaxy_mask_map[unmasked_pixel] = hdu[1].data["FRACGOOD"][index]

    unmasked_pixels = np.where(RM_galaxy_mask_map != 0.0)[0]

    ra_min = 0.0
    ra_max = 2*np.pi
    dec_min = - np.pi/2.0
    dec_max = np.pi/2.0

    ras_of_pixels = np.random.uniform(ra_min, ra_max, len(ras)*60)       
    sindecs_of_pixels = np.random.uniform(np.sin(dec_min), np.sin(dec_max), len(ras)*60)
    decs_of_pixels = np.arcsin(sindecs_of_pixels)

    phis_of_pixels = ras_of_pixels
    thetas_of_pixels = - decs_of_pixels + np.pi/2.0
    pixels = hp.pixelfunc.ang2pix(hp.npix2nside(len(RM_galaxy_mask_map)), thetas_of_pixels, phis_of_pixels)

    proper_rand_ras = ras_of_pixels[np.isin(pixels,unmasked_pixels)]
    proper_rand_decs = decs_of_pixels[np.isin(pixels,unmasked_pixels)]
    rand_galaxy = tr.Catalog(ra=proper_rand_ras, dec=proper_rand_decs, ra_units='radians', dec_units='radians')





    if data_type =="log_realiz":
        photon_map_unmasked = hp.read_map("{0}/../sampled_gamma_ray_map_{1}_{2}_MeV.fits".format(super_core_directory, Emin, Emax))
    else:
        photon_map_unmasked = hp.read_map("{0}/flux_9years_{0}_C.fits".format(gamma_ray_map_unmasked_directory, Emin, Emax))
    if data_type =="log_realiz":
        photon_mask_map_high_res = copy.deepcopy(RM_galaxy_mask_map)
    else:
        photon_mask_map_high_res = hp.read_map("{0}/mask_GP30.0_sources_variable_FL8Y_incl_3FHL_incl_{1}_{2}_MeV_hpx_ord10.fits".format(gamma_ray_mask_map_directory, Emin, Emax))
    photon_map_high_res = photon_map_unmasked * hp.ud_grade(photon_mask_map_high_res, hp.npix2nside(len(photon_map_unmasked)))
    photon_map = hp.ud_grade(photon_map_high_res, nside)
    photon_mask_map = hp.ud_grade(photon_mask_map_high_res, nside)



    nonzero_pixels_photon = np.where(photon_map != 0.0)[0]
    proper_cat_ras_photon, proper_cat_decs_photon, proper_cat_weights_photon = get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_pixels_photon, photon_map)
    cat_photon = tr.Catalog(ra=proper_cat_ras_photon, dec=proper_cat_decs_photon, w=proper_cat_weights_photon, ra_units='radians', dec_units='radians')



    unmasked_pixels_photon = np.where(photon_mask_map != 0.0)[0]
    proper_rand_ras_photon, proper_rand_decs_photon, proper_rand_weights_photon = get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(unmasked_pixels_photon, photon_mask_map)
    rand_photon = tr.Catalog(ra=proper_rand_ras_photon, dec=proper_rand_decs_photon, w=proper_rand_weights_photon, ra_units='radians', dec_units='radians')



    dd_photon_galaxy = tr.NNCorrelation(min_sep=5.0, max_sep=600.0, nbins=12, sep_units='arcminutes')#, bin_slop=0.01)
    dd_photon_galaxy.process(cat_photon, cat_galaxy)
    dr_photon_galaxy = tr.NNCorrelation(min_sep=5.0, max_sep=600.0, nbins=12, sep_units='arcminutes')#, bin_slop=0.01)
    dr_photon_galaxy.process(cat_photon, rand_galaxy)
    rd_photon_galaxy = tr.NNCorrelation(min_sep=5.0, max_sep=600.0, nbins=12, sep_units='arcminutes')#, bin_slop=0.01)
    rd_photon_galaxy.process(rand_photon, cat_galaxy)
    rr_photon_galaxy = tr.NNCorrelation(min_sep=5.0, max_sep=600.0, nbins=12, sep_units='arcminutes')#, bin_slop=0.01)
    rr_photon_galaxy.process(rand_photon, rand_galaxy)

    xi, varxi = dd_photon_galaxy.calculateXi(rr_photon_galaxy, dr_photon_galaxy, rd_photon_galaxy)
    r = np.exp(dd_photon_galaxy.meanlogr)
    print("r (in arcminutes): ")
    print(r)
    print("")
    print("")
    print("xi: ")
    print(xi)
    print("r (in arcminutes) as a list: ")
    print(r.tolist())
    print("xi as a list: ")
    print(xi.tolist())
    np.save("{0}/xi.npy".format(core_directory), xi)
