
logdir=/nfs/slac/g/ki/ki23/des/troughs/flask/log

for ((i=1; i<=1; i++))
do
  l=`cat $logdir/oliver${i}.stdout | grep 
  #bsub -R "rusage[mem=20000]" -eo $logdir/oliver${i}.stderr -oo $logdir/oliver${i}.stdout -W 47:00 bash /u/ki/dgruen/work/flaskpipe/run.sh $i
done

