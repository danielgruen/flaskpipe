#!/bin/bash
# one script to rule them all
source ~/.bash_profile
set +o noglob

#####
# (1) configuration

# (1a) catalogs and input masks to be present in data/

# tracer catalog
tracers=$WORK/fp_${RUN}_$1/${RUN}/catalog.fits
# FITS column names for theta, phi coordinates
tracers_thetacol="dec" # yes, this is a confusing convention in flask
tracers_phicol="ra"

mask_nside=4096
ncookie=$NCOOKIES # number of cookies cut

# (1b) setting for trough catalog
 
trough_nside=1024

# (1d) settings for this script
base=$WORK/fp_${RUN}_$1/
helpers="$PREFIX/trough_helpers/"


##########################
 echo "(1.2) convert tracer catalog to txt file"
   if [ ! -s $base/tracers.txt ]
   then
     echo $helpers/cut_mask $tracers $tracers_phicol $tracers_thetacol $PREFIX/$MASK
     $helpers/cut_mask $tracers $tracers_phicol $tracers_thetacol $PREFIX/$MASK > $base/tracers.txt
   fi
 
##########################

##########################
 echo "(2.1) generate count healpix map from tracer catalog"
 for theta in $countradii
 do
    if [ ! -f $base/trough_${theta}_count.fits* ]
    then
      echo $helpers/count_in_trough $base/tracers.txt $theta $base/trough_${theta}_count.fits
      $helpers/count_in_trough $base/tracers.txt $theta $base/trough_${theta}_count.fits
    fi
 done
##########################



##########################
 echo "(2.3) generate trough maps"
 for(( i=1; i<=$ncookie; ++i ))
 do
  for theta in $countradii
  do
    if [ ! -f $base/trough_${theta}_${i}_0.fits* ]
    then
      echo $helpers/troughfinder $base/trough_${theta}_count.fits* $PREFIX/${COOKIEMASKPREFIX}_c${i}_${theta}.fits.gz ${count_maxmaskfrac} $base/trough_${theta}_$i $countpercentiles
      $helpers/troughfinder $base/trough_${theta}_count.fits* $PREFIX/${COOKIEMASKPREFIX}_c${i}_${theta}.fits.gz ${count_maxmaskfrac} \
                            $base/trough_${theta}_$i $countpercentiles &> $base/trough_${theta}_${i}.log
    fi
  done
 done
##########################

##########################
 echo "(2.4) calculate p(N) statistics"

 for(( i=1; i<=$ncookie; ++i ))
 do
  for theta in $countradii
  do
    if [ ! -s $base/pofn_${theta}_$i.tab ]
    then
      echo $helpers/pofn_bernoulli $base/trough_${theta}_count.fits* ${PREFIX}/${COOKIEMASKPREFIX}_c${i}_${theta}.fits.gz ${count_maxmaskfrac} 
      $helpers/pofn_bernoulli $base/trough_${theta}_count.fits* ${PREFIX}/${COOKIEMASKPREFIX}_c${i}_${theta}.fits.gz \
                             ${count_maxmaskfrac}  > $base/pofn_${theta}_${i}.tab
    fi
  done
 done

##########################

##########################
 echo "(2.4) measure shear"

 for m in f2z1 f2z2
 do
  for theta in $countradii
  do
   for(( i=1; i<=$ncookie; ++i ))
   do
   for(( j=1; i<=$NSOURCEBINS; ++i ))
   do
   for s in $percentilekeys
   do
      if [ ! -s $base/gammat_${m}_${theta}_${s}_${i}_${j}.tab ]
      then
        echo python ${PREFIX}/measure_gammat.py $base/trough_${theta}_${i}_${s}.fits* $base/$RUN/kappa-gamma-${m}.fits* $PREFIX/${COOKIEMASKPREFIX}_c${i}_s${j}.fits.gz $base/gammat_${m}_${theta}_${s}_${i}_$j.tab
        python ${PREFIX}/measure_gammat.py $base/trough_${theta}_${i}_${s}.fits* $base/$RUN/kappa-gamma-${m}.fits* $PREFIX/${COOKIEMASKPREFIX}_c${i}_s${j}.fits.gz $base/gammat_${m}_${theta}_${s}_${i}_$j.tab &> $base/gammat_${m}_${theta}_${s}_${i}_${j}.log
      fi
      #if [ ! -s $base/gammat_bigmask_${m}_${theta}_${s}_${i}_${j}.tab ]
      #then
      #  echo python ${PREFIX}/measure_gammat.py $base/trough_${theta}_${i}_${s}.fits* $base/$RUN/kappa-gamma-${m}.fits* $PREFIX/${COOKIEMASKPREFIX}_c${i}.fits.gz $base/gammat_bigmask_${m}_${theta}_${s}_${i}_$j.tab 
      #  python ${PREFIX}/measure_gammat.py $base/trough_${theta}_${i}_${s}.fits* $base/$RUN/kappa-gamma-${m}.fits* $PREFIX/${COOKIEMASKPREFIX}_c${i}.fits.gz $base/gammat_bigmask_${m}_${theta}_${s}_${i}_$j.tab &> $base/gammat_bigmask_${m}_${theta}_${s}_${i}_${j}.log
      #fi

    done
   done
   done
  done
 done

##########################

##########################
# echo "(2.5) measure shape noise shear"

#  for theta in $countradii
#  do
#   for(( i=1; i<=$ncookie; ++i ))
#   do
#   for s in $percentilekeys
#   do
#      if [ ! -s $base/gammat_shapenoise_${theta}_${s}_$i.tab ]
#      then
#        echo python ${PREFIX}/measure_gammat_shapenoise.py $base/trough_${theta}_${i}_${s}.fits* 0.037 ${PREFIX}/${COOKIEMASKPREFIX}_c${i}_sFIXMEEEEE.fits.gz $base/gammat_shapenoise_${theta}_${s}_$i.tab ${1}$i
#        python ${PREFIX}/measure_gammat_shapenoise.py $base/trough_${theta}_${i}_${s}.fits* 0.037 ${PREFIX}/$COOKIEMASKPREFIX${i}.fits.gz $base/gammat_shapenoise_${theta}_${s}_$i.tab ${1}0$i
#      fi
#   done
#   done
#  done
##########################
