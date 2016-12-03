#!/bin/bash

source /u/ki/dgruen/work/trough_flask/config.rc

if [ $# -ne 1 ]
then
 echo "syntax: $0 ID" 
 exit 1
fi

# (1) make scratch directory
echo mkdir -p $WORK/tf$1
mkdir -p $WORK/tf$1

# (2) copy config file, apply bias
echo cp $PREFIX/${RUN}.config $WORK/tf$1/
cp $PREFIX/${RUN}.config $WORK/tf$1/
echo $PREFIX/flask_bias_info.py ${RUN}-info.dat $BIAS $WORK/tf$1/${RUN}-info.dat
$PREFIX/flask_bias_info.py ${RUN}-info.dat $BIAS $WORK/tf$1/${RUN}-info.dat

# (3) copy Cl files, apply bias
echo $PREFIX/flask_bias_cl.py ${RUN}-info.dat ${RUN}Cl- $BIAS $WORK/tf$1/${RUN}Cl-
$PREFIX/flask_bias_cl.py ${RUN}-info.dat ${RUN}Cl- $BIAS $WORK/tf$1/${RUN}Cl-

# (4) copy n(z) files
echo cp $PREFIX/${RUN}pz-f?.dat  $WORK/tf$1/
cp $PREFIX/${RUN}pz-f?.dat  $WORK/tf$1/

# (5) copy survey mask
echo cp $PREFIX/$MASK $WORK/tf$1/
cp $PREFIX/$MASK $WORK/tf$1/
echo gunzip $WORK/tf$1/$MASK
gunzip $WORK/tf$1/$MASK
echo mv $WORK/tf$1/`basename $MASK .gz` $WORK/tf$1/${RUN}mask-f1z1.fits
mv $WORK/tf$1/`basename $MASK .gz` $WORK/tf$1/${RUN}mask-f1z1.fits


