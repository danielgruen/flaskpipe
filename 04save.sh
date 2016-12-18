#!/bin/bash


if [ $# -ne 2 ]
then
 echo "syntax: $0 ID FLASKPIPE_DIR" 
 exit 1
fi

source $2/config.rc


src=$WORK/fp_${RUN}_$1/$RUN

dest=$SAVEDIR/${RUN}_$1
mkdir -p $dest

# save Poisson-sampled tracer catalogs
cp $src/catalog.fits $dest

# save lensing maps after masking / zipping
for i in $src/kappa-gamma-f?z?.fits
do
python $PREFIX/mask_maps.py $i $PREFIX/$MASK 3 $dest/`basename $i .fits`.fits.gz
done
for i in $src/map-f?z?.fits
do
python $PREFIX/mask_maps.py $i $PREFIX/$MASK 1 $dest/`basename $i .fits`.fits.gz
done

# zip whichever fits might not be zipped yet
for i in $dest/*.fits
do
  gzip $i
done

# clean up scratch directory
rm -rf $WORK/fp_${RUN}_$1/
