#!/bin/bash


if [ $# -ne 2 ]
then
 echo "syntax: $0 ID FLASKPIPE_DIR" 
 exit 1
fi
source $2/config.rc
source $2/config_${RUN}.rc

echo "01CPINPUTS..............."

# (1) make scratch directory
scratchdir=$WORK/fp_${RUN}_$1
echo mkdir -p $scratchdir
mkdir -p $scratchdir

# (2) copy run_flask_and_unscale_maps.py into scratch directory
scratchdir=$WORK/fp_${RUN}_$1
echo cp $PREFIX/run_flask_and_unscale_maps.py $scratchdir
cp $PREFIX/run_flask_and_unscale_maps.py $scratchdir

# (3) copy config file, apply bias
echo cp $PREFIX/${RUN}.config $scratchdir
cp $PREFIX/${RUN}.config $scratchdir
echo $PYTHON $PREFIX/flask_bias_info.py ${RUN}-info.dat ${BIAS[*]} $scratchdir/${RUN}-info.dat
$PYTHON $PREFIX/flask_bias_info.py ${RUN}-info.dat ${BIAS[*]} $scratchdir/${RUN}-info.dat

# (4) copy Cl files, apply bias
echo $PYTHON $PREFIX/flask_bias_cl.py ${RUN}-info.dat ${RUN}Cl- ${BIAS[*]} $scratchdir/${RUN}Cl-
$PYTHON $PREFIX/flask_bias_cl.py ${RUN}-info.dat ${RUN}Cl- ${BIAS[*]} $scratchdir/${RUN}Cl-

# (5) copy n(z) files
#echo cp $PREFIX/${RUN}pz-f?.dat  $scratchdir
#cp $PREFIX/${RUN}pz-f?.dat  $scratchdir

# (6) copy survey mask; this assumes it's the same mask for all fields
echo cp $PREFIX/$MASK $scratchdir
cp $PREFIX/$MASK $scratchdir
echo gunzip $scratchdir/$MASK
gunzip $scratchdir/$MASK
echo mv $scratchdir/`basename $MASK .gz` $scratchdir/mask_to_copy.fits
mv $scratchdir/`basename $MASK .gz` $scratchdir/mask_to_copy.fits
# this does not allow for more than one trough redshift range, but hey

# (7) make copies of mask for all fields
echo mkdir $scratchdir/temp
mkdir $scratchdir/temp
echo cp $scratchdir/mask_to_copy.fits $scratchdir/temp
cp $scratchdir/mask_to_copy.fits $scratchdir/temp


echo $PYTHON $PREFIX/make_copies_of_mask_for_all_fields.py --info_dat_filename $scratchdir/${RUN}-info.dat --directory_with_mask_to_copy $scratchdir --run_name ${RUN}
$PYTHON $PREFIX/make_copies_of_mask_for_all_fields.py --info_dat_filename $scratchdir/${RUN}-info.dat --directory_with_mask_to_copy $scratchdir --run_name ${RUN}

# (8) scale all fields and save the scale values
echo $PYTHON $PREFIX/scale_fields.py --info_dat_filename $scratchdir/${RUN}-info.dat --directory_with_mask_to_copy $scratchdir --run_name ${RUN}
$PYTHON $PREFIX/scale_fields.py --info_dat_filename $scratchdir/${RUN}-info.dat --directory_with_mask_to_copy $scratchdir --run_name ${RUN}


