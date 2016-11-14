#!/bin/bash
# one script to rule them all

set +o noglob

#####
# (1) configuration

# (1a) catalogs and input masks to be present in data/

# tracer catalog
tracers=$WORK/tf$1/$RUN/catalog.fits
# FITS column names for RA, DEC, z coordinates
tracers_thetacol="THETA"
tracers_phicol="PHI"

mask_nside=4096
ncookie=8 # number of cookies cut

# (1b) setting for trough catalog

# count config
countradii="5 10 20" # radii (arcmin) used
countpercentiles="0.2 0.8"  # lower percentile for trough selection
count_maxmaskfrac="0.2" 
trough_nside=1024

# (1d) settings for this script
base=$WORK/tf$1/$RUN
tf="/home/ki/dgruen/trough_flask/"
helpers="/home/ki/dgruen/trough_flask/trough_helpers/"


##########################
 echo "(1.2) convert tracer catalog to txt file"
 for(( i=0; i<$ncookie; ++i ))
 do
   if [ ! -f $base/tracers_$i.txt ]
   then
     echo $helpers/cut_mask $tracers $tracers_thetacol $tracers_phicol $tf/m$i.fits.gz
     $helpers/cut_mask $tracers $tracers_thetacol $tracers_phicol $tf/m$i.fits.gz > $base/tracers_$i.txt
   fi
 done


exit
 
##########################

##########################
 echo "(1.3) convert tracer catalog to healpix map"
 for(( i=0; i<$ncookie; ++i ))
 do
   if [ ! -f $base/tracers_$i.fits ]
   then
      echo $helpers/count_in_trough $base/tracers_$i.txt 0 $base/tracers_$i.fits
      $helpers/count_in_trough $base/tracers_$i.txt 0 $base/tracers_$i.fits
   fi
 done
 
##########################

##########################
 echo "(2.1) generate count healpix map from tracer catalog"
 for theta in $countradii
 do
  for(( i=0; i<$ncookie; ++i ))
  do
    if [ ! -f $base/trough_${theta}_count_$i.fits ]
    then
      echo $helpers/count_in_trough $base/tracers_$i.txt $theta $base/trough_${theta}_count_$i.fits
      $helpers/count_in_trough $base/tracers_$i.txt $theta $base/trough_${theta}_count_$i.fits
    fi
  done
 done
##########################


##########################
 echo "(2.3) generate trough maps"
 for theta in $countradii
 do
  for(( i=0; i<$ncookie; ++i ))
  do
    if [ ! -f $base/trough_${theta}_${lowz}_${highz}_0_$i.fits ]
    then
      echo $helpers/troughfinder trough_${theta}_count_$i.fits trough_${theta}_mask_$i.fits ${count_maxmaskfrac} $base/trough_${theta}_$i $countpercentiles
      $helpers/troughfinder trough_${theta}_count_$i.fits trough_${theta}_mask_$i.fits ${count_maxmaskfrac} $base/trough_${theta}_$i $countpercentiles > trough_${theta}_$i.log
    fi
  done
 done
##########################




















##########################
 echo "(2.4) calculate p(N) statistics"
 mkdir -p pofn
 mkdir -p pofn_jackknife
 cd pofn

 if [ "$incremental" == "0" ]
 then
     rm -f pofn_*.tab
 fi

 for theta in $countradii
 do
  for(( i=0; i<${#countredshift_min[@]}; ++i ))
  do
    lowz=${countredshift_min[$i]}
    highz=${countredshift_max[$i]}

    if [ ! -s pofn_${theta}_${lowz}_${highz}.tab ]
    then
      $helpers/pofn_bernoulli_jackknife ../lenses/trough_${theta}_${lowz}_${highz}_count.fits ../lenses/trough_${theta}_${highz}_mask.fits \
                                             ../lenses/jackknife.fits ${count_maxmaskfrac} ../pofn_jackknife/pofn_${theta}_${lowz}_${highz} > pofn_${theta}_${lowz}_${highz}.tab
    fi
  done
 done

 cd ..
##########################

##########################
 echo "(2.5) calculate p(N1,N2) statistics"
 cd pofn

 if [ "$incremental" == "0" ]
 then
     rm -f pofn1n2_*.tab
 fi

 for theta1 in $countradii 
 do
   for theta2 in $countradii
   do
     if (( $theta1 >= $theta2 ))
     then
       continue
     fi

     for(( i=0; i<${#countredshift_min[@]}; ++i ))
     do
       lowz=${countredshift_min[$i]}
       highz=${countredshift_max[$i]}

       if [ ! -s pofn1n2_${theta1}_${theta2}_${lowz}_${highz}.tab ]
       then
         $helpers/pofn_bernoulli_jackknife_tworadii ../lenses/trough_${theta1}_${lowz}_${highz}_count.fits ../lenses/trough_${theta1}_${highz}_mask.fits $theta1 \
                                                    ../lenses/trough_${theta2}_${lowz}_${highz}_count.fits ../lenses/trough_${theta2}_${highz}_mask.fits $theta2 \
                                                    ../lenses/jackknife.fits ${count_maxmaskfrac} > pofn1n2_${theta1}_${theta2}_${lowz}_${highz}.tab
       fi
     done
   done
 done

 cd ..
##########################



##########################
 
if [ "1" -eq "2" ]
then
 echo "(2.5) make pretty pictures"

 cd lenses
 for theta in $countradii
 do
  for(( i=0; i<${#countredshift_min[@]}; ++i ))
  do
    lowz=${countredshift_min[$i]}
    highz=${countredshift_max[$i]}

    if [ "$incremental" == "0" ]
    then
      rm -f trough_${theta}_${lowz}_${highz}_*.png
    fi

    if [ ! -f trough_${theta}_${lowz}_${highz}_count.png ]
    then
      python $helpers/healpixview.py trough_${theta}_${lowz}_${highz}_count.fits trough_${theta}_${lowz}_${highz}_troughmask.fits
      mv output.png trough_${theta}_${lowz}_${highz}_count.png
    fi

    for f in trough_${theta}_${lowz}_${highz}_[0-9].fits
    do
      if [ ! -f `basename $f .fits`.png ]
      then
        echo python $helpers/healpixview.py $i trough_${theta}_${lowz}_${highz}_troughmask.fits
        python $helpers/healpixview.py $f trough_${theta}_${lowz}_${highz}_troughmask.fits
        mv output.png `basename $f .fits`.png
      fi
    done
  done
 done

 cd ..
fi

##########################

fi



if [ "$do_shear" == "1" ]
then
 mkdir -p shear_jackknife
 mkdir -p shear
 cd shear_jackknife
 for theta in $countradii
 do
  for(( i=0; i<${#countredshift_min[@]}; ++i ))
  do
    lowz=${countredshift_min[$i]}
    highz=${countredshift_max[$i]}

    for (( j=0; j<${#shapecat_zbins_min[@]}; j++ ))
    do
      lowzs=${shapecat_zbins_min[$j]}
      highzs=${shapecat_zbins_max[$j]}

      if [ "$incremental" == "0" ]
      then
        rm -f trough_${theta}_${lowz}_${highz}_[0-9]*_${lowzs}_${highzs}_*.tab
        rm -f ../shear/trough_${theta}_${lowz}_${highz}_[0-9]*_${lowzs}_${highzs}.tab
      fi

      for f in ../lenses/trough_${theta}_${lowz}_${highz}_[0-9]*.fits
      do
        nused=0
        for(( p=0; p<$npatch; p++ ))
        do
          if [ ! -f `basename $f .fits`_${lowzs}_${highzs}_${p}.tab ]
          then
            echo python $helpers/measureshear.py $f ../lenses/jackknife.fits $p ${shapecat_base}/shapes_${lowzs}_${highzs}.fits `basename $f .fits`_${lowzs}_${highzs}_${p}.tab
            python $helpers/measureshear.py $f ../lenses/jackknife.fits $p ${shapecat_base}/shapes_${lowzs}_${highzs}.fits `basename $f .fits`_${lowzs}_${highzs}_${p}.tab &
            nused=`expr $nused + 1`
            if [ $nused -eq $ncores ]
            then
              wait
              nused=0
            fi
          fi
        done
        wait
        
        if [ ! -f ../shear/`basename $f .fits`_${lowzs}_${highzs}.tab ]
        then
          python $helpers/combine_patches.py `basename $f .fits`_${lowzs}_${highzs} $npatch 24 > ../shear/`basename $f .fits`_${lowzs}_${highzs}.tab
        fi

      done

    done
    
    

  done
 done


 cd ..
fi






if [ "$do_spice" == "1" ]
then

  mkdir -p spice
  cd spice
  # resample

  for(( m=0; m<${#spice_maps[@]}; ++m ))
  do

    spice_map=${spice_maps[$m]}
    spice_mask=${spice_masks[$m]}
    echo $spice_map $spice_mask
    if [ "$spice_mask" == "-" ] # generate mask automatically
    then
      $healpix_automask $base/data/$spice_map $base/data/automask.fits
      spice_mask="automask.fits"
    fi

    rm -f map.fits mask.fits

  
    for theta in $countradii
    do
     for(( i=0; i<${#countredshift_min[@]}; ++i ))
     do
      lowz=${countredshift_min[$i]}
      highz=${countredshift_max[$i]}
      for(( b=0; b<=4; b++ ))
      do

       if [ -s trough_${theta}_${lowz}_${highz}_${b}_X_`basename ${spice_map} .fits`.cor ]
       then
        continue
       fi

       if [ ! -f map.fits ]
       then
         echo healpix_resample.py $base/data/${spice_map} ${trough_nside} map.fits
         healpix_resample.py $base/data/${spice_map} ${trough_nside} map.fits
       
         echo healpix_resample.py $base/data/${spice_mask} ${trough_nside} mask.fits 
         healpix_resample.py $base/data/${spice_mask} ${trough_nside} mask.fits 
       fi

       echo spice -mapfile $base/lenses/trough_${theta}_${lowz}_${highz}_$b.fits -maskfile $base/lenses/trough_${theta}_${lowz}_${highz}_troughmask.fits -mapfile2 map.fits -maskfile2 mask.fits -corfile trough_${theta}_${lowz}_${highz}_${b}_X_`basename ${spice_map} .fits`.cor
       spice -mapfile $base/lenses/trough_${theta}_${lowz}_${highz}_$b.fits -maskfile $base/lenses/trough_${theta}_${lowz}_${highz}_troughmask.fits -mapfile2 map.fits -maskfile2 mask.fits -corfile trough_${theta}_${lowz}_${highz}_${b}_X_`basename ${spice_map} .fits`.cor
      done
     done
    done

    rm -f $base/data/automask.fits

  done

  cd ..

fi

if [ "$do_systematics" == "1" ]
then

  mkdir -p systematics
  cd systematics
  # resample

  for(( m=0; m<${#systematics_maps[@]}; ++m ))
  do

    spice_map=${systematics_maps[$m]}
    $healpix_automask $base/data/$spice_map $base/data/automask.fits
    spice_mask="automask.fits"

    rm -f map.fits mask.fits
  
    for theta in $countradii
    do
     for(( i=0; i<${#countredshift_min[@]}; ++i ))
     do
      lowz=${countredshift_min[$i]}
      highz=${countredshift_max[$i]}
      for(( b=0; b<=4; b++ ))
      do

       if [ -s trough_${theta}_${lowz}_${highz}_${b}_X_`basename ${spice_map} .fits`.cor ]
       then
        continue
       fi

       if [ ! -f map.fits ]
       then
         echo healpix_resample.py $base/data/${spice_map} ${trough_nside} map.fits
         healpix_resample.py $base/data/${spice_map} ${trough_nside} map.fits
       
         echo healpix_resample.py $base/data/${spice_mask} ${trough_nside} mask.fits 
         healpix_resample.py $base/data/${spice_mask} ${trough_nside} mask.fits 
       fi

       echo spice -mapfile $base/lenses/trough_${theta}_${lowz}_${highz}_$b.fits -maskfile $base/lenses/trough_${theta}_${lowz}_${highz}_troughmask.fits -mapfile2 map.fits -maskfile2 mask.fits -corfile trough_${theta}_${lowz}_${highz}_${b}_X_`basename ${spice_map} .fits`.cor
       spice -mapfile $base/lenses/trough_${theta}_${lowz}_${highz}_$b.fits -maskfile $base/lenses/trough_${theta}_${lowz}_${highz}_troughmask.fits -mapfile2 map.fits -maskfile2 mask.fits -corfile trough_${theta}_${lowz}_${highz}_${b}_X_`basename ${spice_map} .fits`.cor
      done
     done
    done

    rm -f $base/data/automask.fits

  done

  cd ..

fi

