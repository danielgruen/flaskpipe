import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--RNDSEED', default = "12345")
parser.add_argument('--SELEC_SCALE', default = "0.04377") #SDSS, redMaGiC, z=0.2..0.45, 10' and 20' doesn't matter, corrected for 10 per-cent masking
parser.add_argument('--NSIDE', default = "1024")
parser.add_argument('--SHEAR_LMAX', default = "1024")
parser.add_argument('--LRANGE', default = "10 1032")
parser.add_argument('--flask_executable', default = "/nfs/slac/g/ki/ki19/des/aresh/Old_Stuff/Code_Depository/flask_folder/flask/bin/flask")
parser.add_argument('--troughlenser_executable', default = "/nfs/slac/g/ki/ki23/des/dgruen/trough_lenser/trough_flask_input")
parser.add_argument('--run_name', default = "gxg")
parser.add_argument('--countradii', default = "10 20") #radii (arcmin) used
parser.add_argument('--countradius', default = "20") #the one used for kappa_0, theta_0
parser.add_argument('--countpercentiles', default = "0.2 0.4 0.6 0.8") #lower percentile for trough selection
parser.add_argument('--percentilekeys', default = "0 1 2 3 4 troughmask") #suffixes for measuring shear around
parser.add_argument('--count_maxmaskfrac', default = "0.1")
parser.add_argument('--goodfrac', default = "0.9")
options = parser.parse_args()
RNDSEED = options.RNDSEED
SELEC_SCALE = options.SELEC_SCALE
NSIDE = options.NSIDE
SHEAR_LMAX = options.SHEAR_LMAX
LRANGE = options.LRANGE
flask_executable = options.flask_executable
troughlenser_executable = options.troughlenser_executable
run_name = options.run_name
countradii = options.countradii
countradius = options.countradius
countpercentiles = options.countpercentiles
percentilekeys = options.percentilekeys
count_maxmaskfrac = options.count_maxmaskfrac
print("RNDSEED: {0}".format(RNDSEED))
print("SELEC_SCALE: {0}".format(SELEC_SCALE))
print("NSIDE: {0}".format(NSIDE))
print("SHEAR_LMAX: {0}".format(SHEAR_LMAX))
print("LRANGE: {0}".format(LRANGE))
print("flask_executable: {0}".format(flask_executable))
print("troughlenser_executable: {0}".format(troughlenser_executable))
print("run_name: {0}".format(run_name))
print("countradii: {0}".format(countradii))
print("countradius: {0}".format(countradius))
print("countpercentiles: {0}".format(countpercentiles))
print("percentilekeys: {0}".format(percentilekeys))
print("count_maxmaskfrac: {0}".format(count_maxmaskfrac))



def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension



core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()

FIELDS_INFO = "{0}-info.dat".format(run_name)
CL_PREFIX = "{0}Cl-".format(run_name)
MAPFITS_PREFIX = "{0}/map-".format(run_name)
SHEAR_FITS_PREFIX = "{0}/kappa-gamma-".format(run_name)
flask_config_file = "{0}.config".format(run_name)


#make the out_directory
out_directory = core_directory + "/output/{0}_{1}".format(run_name, RNDSEED)
os_command = "mkdir {0}".format(out_directory)
print(os_command)
os.system(os_command)

#copy everything to out_directory
os_command = "cp {0}* {1} {2} scale_fields.py run_flask_and_unscale_maps.py {3}".format(CL_PREFIX, flask_config_file, FIELDS_INFO, out_directory)
print(os_command)
os.system(os_command)

#make directory for saving the output maps in the out_directory
#os_command = "mkdir {0}/{1}".format(out_directory, run_name)
#print(os_command)
#os.system(os_command)

#make dummy info dat file
os_command = "mv {0}/{1} {0}/{1}_dummy".format(out_directory, FIELDS_INFO)
print(os_command)
os.system(os_command)

#clean the copied info dat file
os_command = "grep -v '^#' {0}/{1}_dummy > {0}/{1}_dummy2".format(out_directory, FIELDS_INFO)
print(os_command)
os.system(os_command)
os_command = "grep -v '^$' {0}/{1}_dummy2 > {0}/{1}".format(out_directory, FIELDS_INFO)
print(os_command)
os.system(os_command)

#delete dummy info dat file
os_command = "rm {0}/{1}_dumm*".format(out_directory, FIELDS_INFO)
print(os_command)
os.system(os_command)

#launch next script
os_command = "python {0}/scale_fields.py --FIELDS_INFO {1} --run_name {2} --flask_executable {3} --RNDSEED {4} --SELEC_SCALE {5} --NSIDE {6} --LRANGE '{7}'".format(out_directory, FIELDS_INFO, run_name, flask_executable, RNDSEED, SELEC_SCALE, NSIDE, LRANGE)
print(os_command)
os.system(os_command)
