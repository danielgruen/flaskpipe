#!/bin/bash


if [ $# -ne 2 ]
then
 echo "syntax: $0 ID FLASKPIPE_DIR" 
 exit 1
fi
source $2/config.rc
source $2/config_${RUN}.rc

echo "02RUNFLASK..............."

cd $WORK/fp_${RUN}_$1
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
       LRANGE: 2 $LMAX \
       MAP_OUT: 0 \
       MAPFITS_PREFIX: $RUN/matter- \
       RECOVCLS_OUT: 0 \
       SHEAR_FITS_PREFIX: $RUN/kappa-gamma- \
       MAPWERFITS_PREFIX: 0 \
       CATALOG_OUT: $RUN/catalog.fits  

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
       LRANGE: 2 $LMAX \
       MAP_OUT: 0 \
       MAPFITS_PREFIX: $RUN/map- \
       RECOVCLS_OUT: 0 \
       SHEAR_FITS_PREFIX: $RUN/kappa-gamma- \
       MAPWERFITS_PREFIX: 0 \
       CATALOG_OUT: $RUN/catalog.fits

