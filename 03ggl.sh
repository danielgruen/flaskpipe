#!/bin/bash


if [ $# -ne 2 ]
then
 echo "syntax: $0 ID FLASKPIPE_DIR" 
 exit 1
fi

source $2/config.rc

echo "03GGL................"

cd $WORK/fp_${RUN}_$1/$RUN

source $PREFIX/ggl_helpers/do_it_all.sh 

# save
src=$WORK/fp_${RUN}_$1/$RUN

dest=$SAVEDIR/${RUN}_$1
mkdir -p $dest

echo "COPYING RESULTS TO $dest ......................"
cp $src/gammat*.tab $dest
cp $src/../gammat*.tab $dest
echo "DONE"
