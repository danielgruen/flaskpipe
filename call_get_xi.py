import os
import glob


def find_directory_for_this_file_and_the_name_of_this_file_without_the_extension():
    file_path = os.path.realpath(__file__)
    program_name = file_path.split("/")[-1]
    core_directory = file_path.split("/{0}".format(program_name))[0]
    program_name_no_extension = program_name.split(".")[0]
    return core_directory, program_name_no_extension



core_directory, program_name_no_extension = find_directory_for_this_file_and_the_name_of_this_file_without_the_extension()



#launch the jobs
directories = glob.glob(core_directory + "/output/*")
for out_directory in directories:
    map_filenames = glob.glob("{0}/map-f*".format(out_directory))
    map_filenames = sorted(map_filenames)
    #copy the scripts and launch the jobs
    for i, map_filename_i in enumerate(map_filenames):
        for j, map_filename_j in enumerate(map_filenames):
            if i>j:
                continue

            map_filename_i_end = map_filename_i.split("map-")[1].split(".")[0]
            map_filename_j_end = map_filename_j.split("map-")[1].split(".")[0]
            #if map_filename_i_end[1] != "6" or map_filename_j_end[1] != "6":
            #    continue
            os_command = "bsub -W 2880 -R bubble -o {0}/get_xi_{1}_{2}.txt python {3}/get_xi.py --map_filename_end_i {1} --map_filename_end_j {2} --out_directory {0}".format(out_directory, map_filename_i_end, map_filename_j_end, core_directory)
            print(os_command)
            os.system(os_command)
