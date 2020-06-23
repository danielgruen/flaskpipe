#!/bin/bash
# change this to reflect path to flaskpipe
d=/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe

if [ $# -ne 1 ]
then
  echo "syntax: run.sh ID"
  echo "where ID is a run number (used as a random seed and for saving outputs)"
  exit 1
fi

echo --------------------------------- bash $d/01cpinputs.sh $1 $d
bash $d/01cpinputs.sh $1 $d
echo --------------------------------- bash $d/02runflask.sh $1 $d
bash $d/02runflask.sh $1 $d

# optional: add your processing steps here to measure and save your data vectors

echo --------------------------------- bash $d/03save.sh $1 $d
bash $d/03save.sh $1 $d
