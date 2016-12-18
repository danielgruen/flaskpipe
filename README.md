This is a set of scripts for creating many log-normal realizations of correlated matter/shear fields over the full sky with flask (https://github.com/hsxavier/flask).

1. Installation
You need flask https://github.com/hsxavier/flask and a python environment with some common (numpy) and less common (healpy) packages.

2. Configuration
To describe your setup of fields and their mutual cross-correlations,
* edit config.rc to contain the correct directories, 
                            bias/density for your galaxy fields,
                            name of healpix mask (possibly containing several realizations in one sky),
                            prefix to healpix masks of individual realizations (only required for optional measurement scripts),
                            NSIDE of mask/output maps and LMAX of your input Cl files
                            a RUN name for your configuration
* place these masks of individual realizations and their stack in this directory in healpix .fits.gz format
  some tools for creating masks: see below
* edit ${RUN}-info.tab to list all galaxy fields (I do that as 1 1 ..., 1 2 ..., 1 3 ...) and shear fields (as 2 1 ..., 2 2 ..., ...)
  see example-info.tab for an example for two foreground fields and to lensing source bins
* copy example.config to ${RUN}.config; likely you don't have to edit anything here (many options are set by 02runflask.sh on the fly)
* create a ${RUN}Cl-XY.dat file containing the C_l for each combination of fields X(=f?z?) and Y in your setup
  see exampleCl-* for an example; note that LMAX in config.rc and the maximum ell in these files must match
* create a ${RUN}pz-f1.dat - just copy the example if you never care for the actual redshift of individual sources in any of your fields

3. Running
run.sh shows how to produce a single realization of outputs. As is usually a good idea for job queue systems, this will copy over files to a work directory, run flask there, and copy back catalogs/maps of density, convergence and shear for storage.
If you want to make any measurements of data vectors, add these between steps 02runflask.sh and 04save.sh.

4. Creating masks
Starting from a single healpix mask of your survey, first resample it to your desired Nside with
healpix_resample.py
you can generate a few rotated masks with
healpix_freerotate.py
and co-add the rotated masks with
healpix_coadd.py

To generate 8 rotations of a DES Y1 SPT footprint m1.fits in a single sky, e.g. use
healpix_freerotate.py m1.fits 180 0 m2.fits
healpix_freerotate.py m1.fits 0 60 m3.fits
healpix_freerotate.py m1.fits 90 60 m4.fits
healpix_freerotate.py m1.fits 180 60 m5.fits
healpix_freerotate.py m1.fits 270 60 m6.fits
healpix_freerotate.py m1.fits 0 120 m7.fits
healpix_freerotate.py m1.fits 180 120 m8.fits


