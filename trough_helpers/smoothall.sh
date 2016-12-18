for ((i=2; i<=8; i++))
do
for t in 30
do
./trough_helpers/smooth_map m$i.fits.gz $t 1024 m${i}_${t}.fits
gzip m${i}_${t}.fits
done
done
