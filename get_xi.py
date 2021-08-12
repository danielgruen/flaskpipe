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
from galXgam.utils import compute_xi_and_r_from_multiple_treecorr_objects
from scipy import interpolate
import treecorr as tr
import numpy as np
import healpy as hp
import pandas as pd
import copy
import argparse
np.random.seed(12345)
parser = argparse.ArgumentParser()
parser.add_argument('--map_filename_end_i')
parser.add_argument('--map_filename_end_j')
parser.add_argument('--out_directory')
options = parser.parse_args()
map_filename_end_i = options.map_filename_end_i
map_filename_end_j = options.map_filename_end_j
out_directory = options.out_directory
print("map_filename_end_i: {0}".format(map_filename_end_i))
print("map_filename_end_j: {0}".format(map_filename_end_j))
print("out_directory: {0}".format(out_directory))




# Load ras, decs, and weights for first field
cat_ras_i = np.load(out_directory + "/{0}_cat_ras.npy".format(map_filename_end_i))
cat_decs_i = np.load(out_directory + "/{0}_cat_decs.npy".format(map_filename_end_i))
cat_weights_i = np.load(out_directory + "/{0}_cat_weights.npy".format(map_filename_end_i))

rand_ras_i = np.load(out_directory + "/{0}_rand_ras.npy".format(map_filename_end_i))
rand_decs_i = np.load(out_directory + "/{0}_rand_decs.npy".format(map_filename_end_i))
rand_weights_i = np.load(out_directory + "/{0}_rand_weights.npy".format(map_filename_end_i))



# Load ras, decs, and weights for second field
cat_ras_j = np.load(out_directory + "/{0}_cat_ras.npy".format(map_filename_end_j))
cat_decs_j = np.load(out_directory + "/{0}_cat_decs.npy".format(map_filename_end_j))
cat_weights_j = np.load(out_directory + "/{0}_cat_weights.npy".format(map_filename_end_j))

rand_ras_j = np.load(out_directory + "/{0}_rand_ras.npy".format(map_filename_end_j))
rand_decs_j = np.load(out_directory + "/{0}_rand_decs.npy".format(map_filename_end_j))
rand_weights_j = np.load(out_directory + "/{0}_rand_weights.npy".format(map_filename_end_j))



# Make treecorr catalogues
cat_i = tr.Catalog(ra=cat_ras_i, dec=cat_decs_i, w=cat_weights_i, ra_units="degrees", dec_units="degrees")
rand_i = tr.Catalog(ra=rand_ras_i, dec=rand_decs_i, w=rand_weights_i, ra_units="degrees", dec_units="degrees")

cat_j = tr.Catalog(ra=cat_ras_j, dec=cat_decs_j, w=cat_weights_j, ra_units="degrees", dec_units="degrees")
rand_j = tr.Catalog(ra=rand_ras_j, dec=rand_decs_j, w=rand_weights_j, ra_units="degrees", dec_units="degrees")



# Compute xi and save
xi, r = compute_xi_and_r_from_multiple_treecorr_objects(cat_i, rand_i, cat_j, rand_j)
np.save(out_directory + "/xi_{0}_{1}.npy".format(map_filename_end_i, map_filename_end_j), xi)
