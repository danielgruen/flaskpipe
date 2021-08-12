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
nside = int(options.NSIDE)
nside_for_galaxy_mask_map = int(options.nside_for_galaxy_mask_map)
map_filename = options.map_filename
out_directory = options.out_directory
print("nside: {0}".format(nside))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("map_filename: {0}".format(map_filename))
print("out_directory: {0}".format(out_directory))




# Enter Directories Here: ####################
galaxy_catalogue_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine"
galaxy_mask_directory = "/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/redMaGiC_Maps_from_redmine"

# Determine the Bin: ####################
map_filename_end = map_filename.split("map-")[1].split(".")[0]
if map_filename_end == "f1z1":
    zmin = 0.15
    zmax = 0.35
if map_filename_end == "f2z2":
    zmin = 0.35
    zmax = 0.50
if map_filename_end == "f3z3":
    zmin = 0.50
    zmax = 0.65
if map_filename_end == "f4z4":
    zmin = 0.65
    zmax = 0.80
if map_filename_end == "f5z5":
    zmin = 0.80
    zmax = 0.90

# Enter Filenames Here: ####################
galaxy_catalogue_filename = galaxy_catalogue_directory + "/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig.fits"
galaxy_mask_filename = galaxy_mask_directory + "/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_mask.fits"
galaxy_overdensity_map_filename = copy.deepcopy(map_filename)



# Masked Sampling From Flask Maps: ####################
true_ras, true_decs, true_weights = read_galaxy_catalogue(galaxy_catalogue_filename, zmin, zmax)


RM_galaxy_mask_map, number_of_galaxies_per_unmasked_pixel = read_galaxy_mask(galaxy_mask_filename, true_weights, nside, nside_for_galaxy_mask_map)
RM_galaxy_mask_map_low_res = np.around(hp.ud_grade(RM_galaxy_mask_map, nside))


galaxy_overdensity_map = hp.read_map(galaxy_overdensity_map_filename)
sampled_galaxy_map = sample_galaxy_overdensity_map(galaxy_overdensity_map, number_of_galaxies_per_unmasked_pixel, nside)


masked_galaxy_map = sampled_galaxy_map * RM_galaxy_mask_map_low_res
cat_ras, cat_decs, cat_weights = map_to_catalogue_in_equatorial_coordinates(masked_galaxy_map, "equatorial", True)
rand_ras, rand_decs = get_random_catalogue_in_equatorial_coordinates_for_binary_mask("equatorial", RM_galaxy_mask_map, 9806760)


np.save(out_directory + "/{0}_cat_ras.npy".format(map_filename_end), cat_ras)
np.save(out_directory + "/{0}_cat_decs.npy".format(map_filename_end), cat_decs)
np.save(out_directory + "/{0}_cat_weights.npy".format(map_filename_end), cat_weights)

np.save(out_directory + "/{0}_rand_ras.npy".format(map_filename_end), rand_ras)
np.save(out_directory + "/{0}_rand_decs.npy".format(map_filename_end), rand_decs)
np.save(out_directory + "/{0}_rand_weights.npy".format(map_filename_end), np.ones_like(rand_ras))

