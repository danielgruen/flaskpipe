
logdir=/nfs/slac/g/ki/ki23/des/troughs/flask/log

for ((i=50; i<100; i++))
do
  bsub -R "rusage[mem=20000]" -eo $logdir/oliver_inplace_${i}.stderr -oo $logdir/oliver_inplace_${i}.stdout -W 24:00 bash /u/ki/dgruen/work/trough_flask/runinplace.sh $i
done

