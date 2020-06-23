
logdir=/nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/log

for ((i=0; i<500; i++))
do
  bsub -J sdss$i -n 4 -R "span[hosts=1] rusage[mem=16000]"  -eo $logdir/sdss${i}.stderr -oo $logdir/sdss${i}.stdout -W 120 bash /nfs/slac/kipac/fs1/g/des/aresh/Gamma_Ray_x_DES/notebook_scripts_and_outputs/using_raw_data/cross/for_all_z_and_E_bins_new_foreground_subtraction/flask_directory/stefano_cls/flaskpipe/run.sh $i
done

