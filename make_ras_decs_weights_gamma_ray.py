import galXgam
from galXgam.utils import get_rms
from galXgam.utils import gal2equat
from galXgam.utils import catalogue_to_map
from galXgam.utils import map_to_array_of_indices_including_repeats
from galXgam.utils import map_to_catalogue_in_equatorial_coordinates
from galXgam.utils import overdensity_to_mean_counts
from galXgam.utils import sample_galaxy_overdensity_map
from galXgam.utils import overdensity_map_to_density_map
from galXgam.utils import overdensity_map_to_density_map_but_keep_negative_values
from galXgam.utils import overdensity_map_to_density_map_but_keep_negative_values_and_use_rms
from galXgam.utils import read_galaxy_catalogue
from galXgam.utils import read_galaxy_mask
from galXgam.utils import read_and_mask_gamma_ray_map
from galXgam.utils import get_gamma_ray_flux_per_unmasked_pixel
from galXgam.utils import get_gamma_ray_rms_of_unmasked_pixels
from galXgam.utils import compute_treecorr_catalogue_and_pair_counts
from galXgam.utils import compute_xi_and_r_from_treecorr_objects
from galXgam.utils import compute_xi_and_r_from_catalogue_and_random_catalogue
from galXgam.utils import get_random_catalogue_in_equatorial_coordinates_across_an_unmasked_sphere
from galXgam.utils import get_random_catalogue_in_equatorial_coordinates_for_binary_mask
from galXgam.utils import transform_C_ell_to_xi
from scipy import interpolate
import treecorr as tr
import numpy as np
import healpy as hp
import pandas as pd
import copy
import argparse
np.random.seed(12345)
parser = argparse.ArgumentParser()
parser.add_argument('--NSIDE', default = "1024")
parser.add_argument('--nside_for_galaxy_mask_map', default = "4096")
parser.add_argument('--map_filename')
parser.add_argument('--out_directory')
options = parser.parse_args()
NSIDE = int(options.NSIDE)
nside_for_galaxy_mask_map = int(options.nside_for_galaxy_mask_map)
map_filename = options.map_filename
out_directory = options.out_directory
print("NSIDE: {0}".format(NSIDE))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("map_filename: {0}".format(map_filename))
print("out_directory: {0}".format(out_directory))




print("0")

# Enter Directories Here: ########################
gamma_ray_map_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/new_stuff/flux_foresub"
gamma_ray_mask_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/new_stuff/masks"
raw_flux_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/new_stuff/flux_wfore"
foreground_counts_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/new_stuff/foreground_counts"
foreground_flux_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/new_stuff/foreground_flux"
exposure_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/raw_photon_exposure_and_foreground_data/new_stuff/exposures"

print("1")

# Determine the Bin: ####################
map_filename_end = map_filename.split("map-")[1].split(".")[0]
if map_filename_end == "f6z1":
    Emin = 631.0
    Emax = 1202.3
if map_filename_end == "f6z2":
    Emin = 1202.3
    Emax = 2290.9
if map_filename_end == "f6z3":
    Emin = 2290.9
    Emax = 4786.3
if map_filename_end == "f6z4":
    Emin = 4786.3
    Emax = 9120.1
if map_filename_end == "f6z5":
    Emin = 9120.1
    Emax = 17378.0
if map_filename_end == "f6z6":
    Emin = 17378.0
    Emax = 36307.8
if map_filename_end == "f6z7":
    Emin = 36307.8
    Emax = 69183.1
if map_filename_end == "f6z8":
    Emin = 69183.1
    Emax = 131825.7
if map_filename_end == "f6z9":
    Emin = 131825.7
    Emax = 1000000.0

print("2")

# Enter Filenames Here: ####################
gamma_ray_map_filename = gamma_ray_map_directory + "/flux_foresub_12years_{0}_{1}_MeV.fits".format(Emin, Emax)
gamma_ray_mask_filename = gamma_ray_mask_directory + "/masks_12years_{0}_{1}_MeV.fits".format(Emin, Emax)
raw_flux_filename = raw_flux_directory + "/flux_wfore_12years_{0}_{1}_MeV.fits".format(Emin, Emax)
foreground_flux_filename = foreground_flux_directory + "/foreground_flux_12years_{0}_{1}_MeV.fits".format(Emin, Emax)
foreground_counts_filename = foreground_counts_directory + "/foreground_counts_12years_{0}_{1}_MeV.fits".format(Emin, Emax)
exposure_filename = exposure_directory + "/exposures_12years_{0}_{1}_MeV.fits".format(Emin, Emax)
gamma_ray_overdensity_map_filename = copy.deepcopy(map_filename)



print("3")

# Masked Sampling from Flask Maps: ####################
# this factor needs to be multiplied with the exposure map to convert between flux and photons properly
steradians_per_pixel = (4.0*np.pi)/hp.nside2npix(1024)


print("4")

# get the gamma ray mask map
gamma_ray_mask_map = hp.read_map(gamma_ray_mask_filename)


print("5")

# get the raw flux map
raw_flux_map = hp.read_map(raw_flux_filename)


print("6")

# get the foreground flux map
foreground_flux_map = hp.read_map(foreground_flux_filename)


print("7")

# get the foreground counts map
foreground_counts_map = hp.read_map(foreground_counts_filename)


print("8")

# get the exposure map
exposure_map = hp.read_map(exposure_filename)


print("9")

# get the gamma ray overdensity map
gamma_ray_overdensity_map = hp.read_map(gamma_ray_overdensity_map_filename)


print("10")

# get the gamma ray flux per unmasked pixel
raw_flux_per_unmasked_pixel = get_gamma_ray_flux_per_unmasked_pixel(raw_flux_map * gamma_ray_mask_map, gamma_ray_mask_map)
foreground_flux_per_unmasked_pixel = get_gamma_ray_flux_per_unmasked_pixel(foreground_flux_map * gamma_ray_mask_map, gamma_ray_mask_map)
gamma_ray_flux_per_unmasked_pixel = raw_flux_per_unmasked_pixel - foreground_flux_per_unmasked_pixel


print("11")

# overdensity to density gamma ray flux map
gamma_ray_density_map = overdensity_map_to_density_map_but_keep_negative_values(gamma_ray_overdensity_map, gamma_ray_flux_per_unmasked_pixel)


print("12")

# turn the density gamma ray flux map to a density gamma ray counts map
photon_map = gamma_ray_density_map * exposure_map * steradians_per_pixel * (Emax - Emin)


print("13")

# add the foreground counts map to the density gamma ray counts map
photon_plus_foreground_map =  photon_map + foreground_counts_map


print("14")

# Poisson sample
photon_plus_foreground_map = np.where(photon_plus_foreground_map >= 0.0, photon_plus_foreground_map, 0.0)
sampled_photon_plus_foreground_map = np.random.poisson(photon_plus_foreground_map).astype(np.float)


print("15")

# convert the Poisson-Sampled Simulated Gamma Ray counts map to a Poisson-Sampled Simulated Gamma Ray flux map
sampled_gamma_ray_plus_foreground_map = sampled_photon_plus_foreground_map * (1.0/exposure_map) * (1.0/steradians_per_pixel) * (1.0/(Emax - Emin))


print("16")

# subtract the foreground
sampled_gamma_ray_map = sampled_gamma_ray_plus_foreground_map - foreground_flux_map


print("17")

# mask the Poisson-Sampled, Foreground-Subtracted Simulated Gamma Ray flux map
sampled_masked_gamma_ray_map = sampled_gamma_ray_map * gamma_ray_mask_map



print("18")

# Get catalogue from masked Poisson-Sampled, Foreground-Subtracted Simulated Gamma Ray flux map: ####################
cat_ras, cat_decs, cat_weights = map_to_catalogue_in_equatorial_coordinates(sampled_masked_gamma_ray_map, "galactic", False)



print("19")

# Using Treecorr to get CF from Simulated Gamma Ray Map
rand_ras, rand_decs = get_random_catalogue_in_equatorial_coordinates_for_binary_mask("galactic", gamma_ray_mask_map, 9806760)


print("20")

# Sace the ras, decs, and weights
np.save(out_directory + "/{0}_cat_ras.npy".format(map_filename_end), cat_ras)
np.save(out_directory + "/{0}_cat_decs.npy".format(map_filename_end), cat_decs)
np.save(out_directory + "/{0}_cat_weights.npy".format(map_filename_end), cat_weights)

np.save(out_directory + "/{0}_rand_ras.npy".format(map_filename_end), rand_ras)
np.save(out_directory + "/{0}_rand_decs.npy".format(map_filename_end), rand_decs)
np.save(out_directory + "/{0}_rand_weights.npy".format(map_filename_end), np.ones_like(rand_ras))

