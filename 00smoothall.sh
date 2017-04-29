#!/bin/bash

if [ $# -ne 1 ]
then
 echo "syntax: $0 FLASKPIPE_DIR" 
 exit 1
fi
source $1/config.rc
source $1/config_${RUN}.rc

for ((i=1; i<=${NCOOKIES}; i++))
do
for t in $countradii
do
if [ ! -f ${COOKIEMASKPREFIX}_c${i}_${t}.fits.gz ]
then
echo ./trough_helpers/smooth_map ${COOKIEMASKPREFIX}$i.fits.gz $t 1024 ${COOKIEMASKPREFIX}${i}_${t}.fits.gz
./trough_helpers/smooth_map ${COOKIEMASKPREFIX}$i.fits.gz $t 1024 ${COOKIEMASKPREFIX}${i}_${t}.fits.gz
fi
done
done
