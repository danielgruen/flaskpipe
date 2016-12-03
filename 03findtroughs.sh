#!/bin/bash

source /u/ki/dgruen/work/trough_flask/config.rc

if [ $# -ne 1 ]
then
 echo "syntax: $0 ID" 
 exit 1
fi

cd $WORK/tf$1/$RUN

source $PREFIX/trough_helpers/do_it_all.sh 
