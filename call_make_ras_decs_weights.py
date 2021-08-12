import os
import argparse
import glob
parser = argparse.ArgumentParser()
parser.add_argument('--number_of_RNDSEED', default = "50")
parser.add_argument('--start_of_RNDSEED', default = "12345")
parser.add_argument('--NSIDE', default = "1024")
parser.add_argument('--nside_for_galaxy_mask_map', default = "4096")
parser.add_argument('--run_name', default = "gxg")
options = parser.parse_args()
number_of_RNDSEED = int(options.number_of_RNDSEED)
start_of_RNDSEED = int(options.start_of_RNDSEED)
NSIDE = options.NSIDE
nside_for_galaxy_mask_map = options.nside_for_galaxy_mask_map
run_name = options.run_name
print("number_of_RNDSEED: {0}".format(number_of_RNDSEED))
print("start_of_RNDSEED: {0}".format(start_of_RNDSEED))
print("NSIDE: {0}".format(NSIDE))
print("nside_for_galaxy_mask_map: {0}".format(nside_for_galaxy_mask_map))
print("run_name: {0}".format(run_name))



def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension



core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()


#make the out_directory and work directories
output_directory = core_directory + "/output"

#copy the scripts and launch the jobs
for i in list(range(0,number_of_RNDSEED)):
    RNDSEED = start_of_RNDSEED + i

    out_directory = core_directory + "/output/{0}_{1}".format(run_name, RNDSEED)
    map_filenames = glob.glob("{0}/map-f*".format(out_directory))
    #copy the scripts and launch the jobs
    for map_filename in map_filenames:
        map_filename_end = map_filename.split("map-f")[1][0]
        map_filename_long_end = map_filename.split("map-")[1].split(".")[0]
        if map_filename_end == "6":
            os_command = "bsub -W 50 -R bubble -o {0}/make_ras_decs_weights_{1}.txt python {2}/make_ras_decs_weights_gamma_ray.py --NSIDE {3} --nside_for_galaxy_mask_map {4} --map_filename {5} --out_directory {0}".format(out_directory, map_filename_long_end, core_directory, NSIDE, nside_for_galaxy_mask_map, map_filename)
            print(os_command)
            os.system(os_command)
        else:
            os_command = "bsub -W 50 -R bubble -o {0}/make_ras_decs_weights_{1}.txt python {2}/make_ras_decs_weights_galaxy.py --NSIDE {3} --nside_for_galaxy_mask_map {4} --map_filename {5} --out_directory {0}".format(out_directory, map_filename_long_end, core_directory, NSIDE, nside_for_galaxy_mask_map, map_filename)
            print(os_command)
            os.system(os_command)
