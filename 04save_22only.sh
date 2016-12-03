#!/bin/bash

source /u/ki/dgruen/work/trough_flask/config.rc

if [ $# -ne 1 ]
then
 echo "syntax: $0 ID" 
 exit 1
fi

DESDIR=/nfs/slac/g/ki/ki23/des/troughs

mkdir -p $DESDIR/flask/${RUN}$1

python $PREFIX/mask_maps.py $WORK/tf$1/$RUN/kappa-gamma-f2z2.fits $PREFIX/$MASK $DESDIR/flask/${RUN}$1/kappa-gamma-f2z2.fits.gz


for i in $DESDIR/flask/${RUN}$1/*.fits
do
  gzip $i
done

