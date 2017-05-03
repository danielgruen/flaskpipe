#!/bin/bash

if [ $# -ne 1 ]
then
 echo "syntax: $0 FLASKPIPE_DIR" 
 exit 1
fi
source $1/config.rc
source $1/config_${RUN}.rc

echo "ONLY WORKS FOR $NSOURCEBINS EQUAL 1"

echo $TROUGHLENSER 0.286 0.82 0 $countradius nofz_redmagic_0.2_0.45_true.dat ${RUN}_s1_pz.tab 0.16 0.64 `calc.sh 3.14159\*$countradius*\$countradius*${DENSITY[0]}\*$goodfrac` ${BIAS[0]} 1 0 cl
#$TROUGHLENSER 0.286 0.82 0 $countradius nofz_redmagic_0.2_0.45_true.dat ${RUN}_s1_pz.tab 0.16 0.64 `calc.sh 3.14159\*$countradius*\$countradius*${DENSITY[0]}\*$goodfrac` ${BIAS[0]} 1 0 cl


# generate info file
rm -f ${RUN}-info.dat
echo "# Field number, z bin number, mean, shift, field type, zmin, zmax" >> ${RUN}-info.dat
echo "# Types: 1-galaxies 2-shear" >> ${RUN}-info.dat
echo >> ${RUN}-info.dat
echo "    1    1   0.0000   `head -n 1 FLASK_input_${countradius}arcmin.tab | cut -d \  -f 4`      1   0.2000   0.4500" >> ${RUN}-info.dat
echo "    2    1   0.0000   `head -n 2 FLASK_input_${countradius}arcmin.tab | tail -n 1 | cut -d \  -f 4`      2   0.2000   0.4500" >> ${RUN}-info.dat
echo "    2    2   0.0000   1000.               2   0.0000   2.0000" >> ${RUN}-info.dat



# generate cl file

tail -n +5 FLASK_input_${countradius}arcmin.tab | awk '{print $1,$2}' > ${RUN}Cl-f1z1f1z1.dat # autocorrelation of matter
tail -n +5 FLASK_input_${countradius}arcmin.tab | awk '{print $1,$3}' > ${RUN}Cl-f1z1f2z1.dat # matter X correlated kappa
tail -n +5 FLASK_input_${countradius}arcmin.tab | awk '{print $1,$4}' > ${RUN}Cl-f2z1f2z1.dat # correlated kappa x correlated kappa
tail -n +5 FLASK_input_${countradius}arcmin.tab | awk '{print $1,$5}' > ${RUN}Cl-f2z2f2z2.dat # autocorrelation of uncorrelated kappa

tail -n +5 FLASK_input_${countradius}arcmin.tab | awk '{print $1,0}' > ${RUN}Cl-f1z1f2z2.dat # matter X uncorrelated kappa
tail -n +5 FLASK_input_${countradius}arcmin.tab | awk '{print $1,0}' > ${RUN}Cl-f2z1f2z2.dat # correlated kappa X uncorrelated kappa


