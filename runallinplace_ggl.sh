
logdir=/nfs/slac/g/ki/ki23/des/troughs/flask/log

for ((i=51; i<150; i++))
do
  bsub -R "rusage[mem=20000]" -eo $logdir/oliver_ggl_${i}.stderr -oo $logdir/oliver_ggl_${i}.stdout -W 24:00 bash /u/ki/dgruen/work/trough_flask/runinplace_ggl.sh $i
done

