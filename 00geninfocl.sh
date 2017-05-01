#!/bin/bash

if [ $# -ne 1 ]
then
 echo "syntax: $0 FLASKPIPE_DIR" 
 exit 1
fi
source $1/config.rc
source $1/config_${RUN}.rc

# generate info file
echo "CANNOT CREATE INFO FILE BECAUSE, WELL, I DO NOT KNOW THE LOGNORMAL PARAMETERS"
cp oliver-info.dat ${RUN}-info.dat

# generate cl file
echo "CANNOT CREATE CL FILES BECAUSE, WELL, I DO NOT KNOW THE POWER SPECTRA"

cp oliverCl-f1z1f1z1.dat ${RUN}Cl-f1z1f1z1.dat # autocorrelation of matter

for ((i=2; i<`expr 2 + $NSOURCEBINS`; i++))
do
  cp oliverCl-f1z1f${i}z1.dat ${RUN}Cl-f1z1f${i}z1.dat # correlated kappa x trough matter
  cp oliverCl-f1z1f${i}z2.dat ${RUN}Cl-f1z1f${i}z2.dat # uncorrelated kappa x trough matter (zeros!)
  cp oliverCl-f${i}z1f${i}z1.dat ${RUN}Cl-f${i}z1f${i}z1.dat # autocorrelation of correlated kappa
  cp oliverCl-f${i}z1f${i}z2.dat ${RUN}Cl-f${i}z1f${i}z2.dat # uncorrelated kappa x correlated kappa (zeros!)
  cp oliverCl-f${i}z2f${i}z2.dat ${RUN}Cl-f${i}z2f${i}z2.dat # autocorrelation of uncorrelated kappa

  for ((j=`expr $i + 1`; j<`expr 2 + $NSOURCEBINS`; j++)) # cross-correlation of source bin kappas
  do
    cp oliverCl-f${i}z1f${j}z1.dat ${RUN}Cl-f${i}z1f${j}z1.dat
    cp oliverCl-f${i}z2f${j}z1.dat ${RUN}Cl-f${i}z2f${j}z1.dat
    cp oliverCl-f${i}z1f${j}z2.dat ${RUN}Cl-f${i}z1f${j}z2.dat
    cp oliverCl-f${i}z2f${j}z2.dat ${RUN}Cl-f${i}z2f${j}z2.dat
  done
done
