#!/bin/bash


if [ $# -ne 2 ]
then
 echo "syntax: $0 ID FLASKPIPE_DIR" 
 exit 1
fi

source $2/config.rc

cd $WORK/fp_${RUN}_$1/$RUN

source $PREFIX/trough_helpers/do_it_all.sh 

# save
src=$WORK/fp_${RUN}_$1/$RUN

dest=$SAVEDIR/${RUN}_$1
mkdir -p $dest

cp $src/pofn*.tab $dest
cp $src/gammat*.tab $dest
cp $src/trough*.fits $dest
