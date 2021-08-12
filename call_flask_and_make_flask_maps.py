import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--number_of_RNDSEED', default = "50")
parser.add_argument('--start_of_RNDSEED', default = "12345")
parser.add_argument('--SELEC_SCALE', default = "0.04377") #SDSS, redMaGiC, z=0.2..0.45, 10' and 20' doesn't matter, corrected for 10 per-cent masking
parser.add_argument('--NSIDE', default = "1024")
parser.add_argument('--SHEAR_LMAX', default = "1024")
parser.add_argument('--LRANGE', default = "10 1032")
parser.add_argument('--flask_executable', default = "/sdf/group/kipac/g/des/ki19/aresh/Old_Stuff/Code_Depository/flask_folder/flask/bin/flask")
parser.add_argument('--run_name', default = "gxg")
options = parser.parse_args()
number_of_RNDSEED = int(options.number_of_RNDSEED)
start_of_RNDSEED = int(options.start_of_RNDSEED)
SELEC_SCALE = options.SELEC_SCALE
NSIDE = options.NSIDE
SHEAR_LMAX = options.SHEAR_LMAX
LRANGE = options.LRANGE
flask_executable = options.flask_executable
run_name = options.run_name
print("number_of_RNDSEED: {0}".format(number_of_RNDSEED))
print("start_of_RNDSEED: {0}".format(start_of_RNDSEED))
print("SELEC_SCALE: {0}".format(SELEC_SCALE))
print("NSIDE: {0}".format(NSIDE))
print("SHEAR_LMAX: {0}".format(SHEAR_LMAX))
print("LRANGE: {0}".format(LRANGE))
print("flask_executable: {0}".format(flask_executable))
print("run_name: {0}".format(run_name))



def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension



core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()


#make the output, work, and log directories
output_directory = core_directory + "/output"
os_command = "mkdir {0}".format(output_directory)
print(os_command)
os.system(os_command)
work_directory = core_directory + "/work"
os_command = "mkdir {0}".format(work_directory)
print(os_command)
os.system(os_command)

#launch the jobs
for i in list(range(0,number_of_RNDSEED)):
    RNDSEED = start_of_RNDSEED + i

    FIELDS_INFO = "{0}-info.dat".format(run_name)
    CL_PREFIX = "{0}Cl-".format(run_name)
    MAPFITS_PREFIX = "{0}/map-".format(run_name)
    SHEAR_FITS_PREFIX = "{0}/kappa-gamma-".format(run_name)
    flask_config_file = "{0}.config".format(run_name)

    #make out and fp directories
    out_directory = core_directory + "/output/{0}_{1}".format(run_name, RNDSEED)
    os_command = "mkdir {0}".format(out_directory)
    print(os_command)
    os.system(os_command)
    fp_directory = core_directory + "/work/fp_{0}_{1}".format(run_name, RNDSEED)
    os_command = "mkdir {0}".format(fp_directory)
    print(os_command)
    os.system(os_command)

    #copy everything to out_directory and one thing to the fp_directory
    os_command = "cp {0}* {1} {2} {3}".format(CL_PREFIX, flask_config_file, FIELDS_INFO, out_directory)
    print(os_command)
    os.system(os_command)
    os_command = "cp {0} {1}".format(flask_config_file, fp_directory)
    print(os_command)
    os.system(os_command)

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



    #launch flask
    os_command = "cd {7}; pwd;bsub -W 50 -R bubble -o {7}/flask.txt {0} {1}.config RNDSEED: {2} FIELDS_INFO: {3} CL_PREFIX: {1}Cl- SELEC_SCALE: {4} NSIDE: {5} SHEAR_LMAX: {5} LRANGE: {6} MAPFITS_PREFIX: map- SHEAR_FITS_PREFIX: kappa-gamma-".format(flask_executable, run_name, RNDSEED, FIELDS_INFO, SELEC_SCALE, NSIDE, LRANGE, out_directory)
    print(os_command)
    os.system(os_command)
