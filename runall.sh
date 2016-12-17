
logdir=/nfs/slac/g/ki/ki23/des/troughs/flask/log

for ((i=200; i<300; i++))
do
  bsub -R "rusage[mem=20000]" -eo $logdir/oliver${i}.stderr -oo $logdir/oliver${i}.stdout -W 47:00 bash /u/ki/dgruen/work/trough_flask/run.sh $i
done

