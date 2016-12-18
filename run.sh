#!/bin/bash
# change this to reflect path to flaskpipe
d=/u/ki/dgruen/work/flaskpipe

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
#echo --------------------------------- bash $d/03findtroughs.sh $1 $d
#bash $d/03findtroughs.sh $1 $d


echo --------------------------------- bash $d/04save.sh $1 $d
bash $d/04save.sh $1 $d
