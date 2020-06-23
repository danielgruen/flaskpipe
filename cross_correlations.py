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



def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension_and_the_energy_range():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    energy_range_preliminary = file_path.split("new_foreground_subtraction_")[1]
    energy_range = energy_range_preliminary.split("/")[0]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension, energy_range



def restrict_generic_map_and_generic_mask_map_using_unmasked_pixels_generic_and_galat_range_and_galon_range(generic_map, generic_mask_map, unmasked_pixels_generic, galat_range=(-np.inf,np.inf), galon_range=(-np.inf,np.inf), abs_galat_min=np.inf):
    galat_min = galat_range[0]
    galat_max = galat_range[1]
    galon_min = galon_range[0]
    galon_max = galon_range[1]

    thetas_of_pixels, phis_of_pixels = hp.pixelfunc.pix2ang(nside,unmasked_pixels_generic)

    galons = []
    galats = []
    for theta_gal, phi_gal, unmasked_pixel in zip(thetas_of_pixels, phis_of_pixels, unmasked_pixels_generic):
        galon = np.degrees(phi_gal)
        galat = -np.degrees(theta_gal) + 90.0
        if galat > galat_max or galat < galat_min or galon > galon_max or galon < galon_min or np.abs(galat) < abs_galat_min:
            generic_map[unmasked_pixel] = 0.0
            generic_mask_map[unmasked_pixel] = 0
        else:
            galons.append(galon)
            galats.append(galat)

    plt.figure()
    plt.scatter(galons,galats)
    plt.xlabel("galon")
    plt.ylabel("galat")
    plt.savefig("{0}/{1}_unmasked_galons_vs_galats.png".format(core_directory, program_name_no_extension))

    return generic_map, generic_mask_map



def get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_or_unmasked_pixels_generic, generic_map_or_generic_mask_map):
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



super_core_directory, program_name_no_extension, energy_range = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension_and_the_energy_range()
zmins = [zbin_edge1, zbin_edge2, zbin_edge3, zbin_edge4, zbin_edge5]
zmaxs = [zbin_edge2, zbin_edge3, zbin_edge4, zbin_edge5, zbin_edge6]

for zmin, zmax in zip(zmins,zmaxs):
    core_directory = super_core_directory + "/z_{0}_{1}".format(zmin,zmax)
    RM_galaxy_mask_map = np.zeros(201326592)
    if data_type =="log_realiz":
        filename = "{0}/../../../../mask_sdss_sum.fits.gz".format(core_directory)
        filename_for_mask = copy.deepcopy(filename)
    else:
        filename = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits"
    hdu = fits.open(filename)
    for index in list(range(0,len(hdu[1].data["HPIX"]))):
        unmasked_pixel = hdu[1].data["HPIX"][index]
        RM_galaxy_mask_map[unmasked_pixel] = hdu[1].data["FRACGOOD"][index]







    unmasked_pixels = []
    for index in list(range(0,201326592)):
        if RM_galaxy_mask_map[index] > 0.0:
            unmasked_pixels.append(index)
    unmasked_pixels = np.array(unmasked_pixels)
    print(len(unmasked_pixels))



    ra_min = 0.0
    ra_max = 2*np.pi
    dec_min = - np.pi/2.0
    dec_max = np.pi/2.0
    number_of_randoms = 0
    ras_randoms = []
    decs_randoms = []

    ras_of_pixels = np.random.uniform(ra_min, ra_max, 22154906)       
    sindecs_of_pixels = np.random.uniform(np.sin(dec_min), np.sin(dec_max), 22154906)
    decs_of_pixels = np.arcsin(sindecs_of_pixels)

    phis_of_pixels = ras_of_pixels
    thetas_of_pixels = - decs_of_pixels + np.pi/2.0
    pixels = hp.pixelfunc.ang2pix(4096, thetas_of_pixels, phis_of_pixels)

    proper_rand_ras = ras_of_pixels[np.isin(pixels,unmasked_pixels)]
    proper_rand_decs = decs_of_pixels[np.isin(pixels,unmasked_pixels)]



    rand_galaxy = tr.Catalog(ra=proper_rand_ras, dec=proper_rand_decs, ra_units='radians', dec_units='radians')


    if data_type =="log_realiz":
        ras = np.load("{0}/../../sampled_ras.npy".format(core_directory))
        decs = np.load("{0}/../../sampled_decs.npy".format(core_directory))
        weights = np.ones(len(ras))
    else:
        filename = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits"
        hdu = fits.open(filename)



        ras = hdu[1].data['RA']
        decs = hdu[1].data['DEC']
        weights = hdu[1].data['weight']
        zs = hdu[1].data['ZREDMAGIC']
        delete_list = []
        for i, z in enumerate(zs):
            if z < zmin or z > zmax:
                delete_list.append(i)
                
        ras = np.delete(ras, delete_list)
        decs = np.delete(decs, delete_list)
        weights = np.delete(weights, delete_list)



    cat_galaxy = tr.Catalog(ra=ras, dec=decs, w=weights, ra_units='degrees', dec_units='degrees')








    if data_type =="log_realiz":
        photon_map_unmasked = hp.read_map("{0}/../../sampled_gamma_ray_map.fits".format(core_directory))
    else:
        photon_map_unmasked = hp.read_map("/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/flux_9years_foresub/flux_9years_{0}_C.fits".format(energy_range))
    if data_type =="log_realiz":
        photon_mask_map_high_res = copy.deepcopy(RM_galaxy_mask_map)
    else:
        photon_mask_map_high_res = hp.read_map("/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/masks/mask_GP30.0_sources_variable_FL8Y_incl_3FHL_incl_{0}_hpx_ord10.fits".format(energy_range))
    photon_map_high_res = photon_map_unmasked * photon_mask_map_high_res



    photon_map = hp.ud_grade(photon_map_high_res, nside)
    photon_mask_map = hp.ud_grade(photon_mask_map_high_res, nside)



    unmasked_pixels_photon = get_unmasked_pixels_generic_given_generic_mask_map(photon_mask_map)
    #photon_map, photon_mask_map = restrict_generic_map_and_generic_mask_map_using_unmasked_pixels_generic_and_galat_range_and_galon_range(photon_map, photon_mask_map, unmasked_pixels_photon, abs_galat_min=0.0)
    #unmasked_pixels_photon = get_unmasked_pixels_generic_given_generic_mask_map(photon_mask_map)



    proper_rand_ras_photon, proper_rand_decs_photon, proper_rand_weights_photon = get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(unmasked_pixels_photon, photon_mask_map)
    rand_photon = tr.Catalog(ra=proper_rand_ras_photon, dec=proper_rand_decs_photon, w=proper_rand_weights_photon, ra_units='radians', dec_units='radians')



    nonzero_pixels_photon = get_nonzero_pixels_generic_given_generic_map(photon_map)
    proper_cat_ras_photon, proper_cat_decs_photon, proper_cat_weights_photon = get_ras_and_decs_and_weights_given_nonzero_or_unmasked_pixels_generic_and_generic_map_or_generic_mask_map(nonzero_pixels_photon, photon_map)
    cat_photon = tr.Catalog(ra=proper_cat_ras_photon, dec=proper_cat_decs_photon, w=proper_cat_weights_photon, ra_units='radians', dec_units='radians')



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



    lmin = 50
    lmax = 1000.0
    lambda_min = (180.0 / lmax) * 60.0
    lambda_max = (180.0 / lmin) * 60.0



    plt.figure()
    plt.plot(r, r*xi)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlim(2.4,230.0)
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel(r'$\theta$' + "w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_CF_short_range.png".format(core_directory, program_name_no_extension))

    plt.figure()
    plt.plot(r, r*xi)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel(r'$\theta$' + "w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_CF_long_range.png".format(core_directory, program_name_no_extension))

    xi_manual = dd_weight_final      /     (rr_weight_final*(dd_tot_final/rr_tot_final))
    plt.figure()
    plt.plot(r, r*xi_manual)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel(r'$\theta$' + "w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_dd_over_rr.png".format(core_directory, program_name_no_extension))

    print("")
    print("")
    print("dd_weight_final      /     (rr_weight_final*(dd_tot_final/rr_tot_final)): ")
    print(dd_weight_final      /     (rr_weight_final*(dd_tot_final/rr_tot_final)))

    xi_manual = rr_weight_final
    plt.figure()
    plt.plot(r, r*xi_manual)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel(r'$\theta$' + "w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_rr.png".format(core_directory, program_name_no_extension))

    print("")
    print("")
    print("rr_weight_final: ")
    print(rr_weight_final)



    plt.figure()
    plt.plot(r, xi)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlim(2.4,230.0)
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel("w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_CF_short_range_without_theta.png".format(core_directory, program_name_no_extension))

    plt.figure()
    plt.plot(r, xi)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel("w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_CF_long_range_without_theta.png".format(core_directory, program_name_no_extension))

    xi_manual = dd_weight_final      /     (rr_weight_final*(dd_tot_final/rr_tot_final))
    plt.figure()
    plt.plot(r, xi_manual)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel("w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_dd_over_rr_without_theta.png".format(core_directory, program_name_no_extension))

    xi_manual = rr_weight_final
    plt.figure()
    plt.plot(r, xi_manual_final)
    plt.axvline(x=lambda_min, label="PSF or photon statisitcs-Ackermann")
    plt.axvline(x=lambda_max, label="foreground emission-Ackermann", color="chocolate")
    plt.xscale("log")
    plt.xlabel(r'$\theta$' + "[arcmin]")
    plt.ylabel("w(" + r'$\theta$' + ")")
    plt.legend()
    plt.savefig("{0}/{1}_rr_without_theta.png".format(core_directory, program_name_no_extension))
