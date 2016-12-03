#!/bin/bash

source /u/ki/dgruen/work/trough_flask/config.rc

if [ $# -ne 1 ]
then
 echo "syntax: $0 ID" 
 exit 1
fi

cd $WORK/tf$1
mkdir -p $RUN

echo $FLASK ${RUN}.config  \
       RNDSEED: $1    \
       ELLIP_SIGMA: 0 \
       FIELDS_INFO: ${RUN}-info.dat \
       CL_PREFIX: ${RUN}Cl- \
       SELEC_PREFIX: ${RUN}mask- \
       SELEC_Z_PREFIX: ${RUN}pz- \
       SELEC_TYPE: 0  \
       SELEC_SEPARABLE: 0  \
       SELEC_SCALE: $DENSITY \
       NSIDE: $NSIDE  \
       SHEAR_LMAX: $NSIDE \
       MAP_OUT: 0 \
       MAPFITS_PREFIX: 0 \
       RECOVCLS_OUT: 0 \
       SHEAR_FITS_PREFIX: $RUN/kappa-gamma- \
       MAPWERFITS_PREFIX: 0 \
       CATALOG_OUT: $RUN/catalog.fits &> $RUN/flask.log 

$FLASK ${RUN}.config  \
       RNDSEED: $1    \
       ELLIP_SIGMA: 0 \
       FIELDS_INFO: ${RUN}-info.dat \
       CL_PREFIX: ${RUN}Cl- \
       SELEC_PREFIX: ${RUN}mask- \
       SELEC_Z_PREFIX: ${RUN}pz- \
       SELEC_TYPE: 0  \
       SELEC_SEPARABLE: 0  \
       SELEC_SCALE: $DENSITY \
       NSIDE: $NSIDE  \
       SHEAR_LMAX: $NSIDE \
       MAP_OUT: 0 \
       MAPFITS_PREFIX: 0 \
       RECOVCLS_OUT: 0 \
       SHEAR_FITS_PREFIX: $RUN/kappa-gamma- \
       MAPWERFITS_PREFIX: 0 \
       CATALOG_OUT: $RUN/catalog.fits &> $RUN/flask.log 

