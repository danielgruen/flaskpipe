Note:

You must add your C_ell files to this directory. For the currently-used C_ells copy over the C_ell files from the "C_ells_unitless_and_splined" directory.


Note:

If you don't specify a flask executable, it will default to using mine



From this directory, run:

1. python call_flask_and_make_flask_maps.py --run_name sdss --number_of_RNDSEED NUMBER_OF_RANDOM_SEEDS_YOU_WANT --flask_executable PATH/TO/FLASK/EXECUTABLE/flask

2. This will launch flask as a batch job (once for each RNDSEEd) and create an "output" and "work" folder. Unless you changed "run_name," there should be subfolders in "output" called "sdss_RNDSEED." In those subfolders, all the outputs will be saved. Wait for flask to finish making overdensity maps (there is a log file for each RNDSEED in each subfolder called "flask.txt").

3. python call_make_ras_decs_and_weights.py --run_name sdss --number_of_RNDSEED NUMBER_OF_RANDOM_SEEDS_YOU_WANT

4. For every RNDSEED, this will launch a batch job that will sample the overdensity map made by flask and make npy files for the ras, decs, and weights for the catalogues and randoms. Wait for all jobs to finish.

5. python call_get_xi.py 

6. This launches batch jobs that will produce a CF for each pair of fields for every RNDSEED. For cases where both fields are gamma ray fields, this could take up to 3 hours or so.
