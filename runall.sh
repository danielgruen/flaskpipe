
logdir=/nfs/slac/g/ki/ki23/des/troughs/flask/log

for ((i=0; i<100; i++))
do
  bsub -J sdss$i -n 4 -R "span[hosts=1] rusage[mem=16000]"  -eo $logdir/sdss${i}.stderr -oo $logdir/sdss${i}.stdout -W 120:00 bash /u/ki/dgruen/work/flaskpipe/run.sh $i
done

