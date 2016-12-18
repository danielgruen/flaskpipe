#!/bin/bash


if [ $# -ne 2 ]
then
 echo "syntax: $0 ID FLASKPIPE_DIR" 
 exit 1
fi
source $2/config.rc

# (1) make scratch directory
scratchdir=$WORK/fp_${RUN}_$1
echo mkdir -p $scratchdir
mkdir -p $scratchdir

# (2) copy config file, apply bias
echo cp $PREFIX/${RUN}.config $scratchdir
cp $PREFIX/${RUN}.config $scratchdir
echo $PYTHON $PREFIX/flask_bias_info.py ${RUN}-info.dat ${BIAS[*]} $scratchdir/${RUN}-info.dat
$PYTHON $PREFIX/flask_bias_info.py ${RUN}-info.dat ${BIAS[*]} $scratchdir/${RUN}-info.dat

# (3) copy Cl files, apply bias
echo $PYTHON $PREFIX/flask_bias_cl.py ${RUN}-info.dat ${RUN}Cl- ${BIAS[*]} $scratchdir/${RUN}Cl-
$PYTHON $PREFIX/flask_bias_cl.py ${RUN}-info.dat ${RUN}Cl- ${BIAS[*]} $scratchdir/${RUN}Cl-

# (4) copy n(z) files
echo cp $PREFIX/${RUN}pz-f?.dat  $scratchdir
cp $PREFIX/${RUN}pz-f?.dat  $scratchdir

# (5) copy survey mask; this assumes it's the same mask for all fields
echo cp $PREFIX/$MASK $scratchdir
cp $PREFIX/$MASK $scratchdir
echo gunzip $scratchdir/$MASK
gunzip $scratchdir/$MASK
echo mv $scratchdir/`basename $MASK .gz` $scratchdir/${RUN}mask-f1z1.fits
mv $scratchdir/`basename $MASK .gz` $scratchdir/${RUN}mask-f1z1.fits


