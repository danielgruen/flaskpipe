
d=/u/ki/dgruen/work/trough_flask
echo --------------------------------- bash $d/01cpinputs.sh $1
bash $d/01cpinputs.sh $1
echo --------------------------------- bash $d/02runflask.sh $1
bash $d/02runflask.sh $1
echo --------------------------------- bash $d/03findtroughs.sh $1
bash $d/03findtroughs.sh $1
echo --------------------------------- bash $d/04save.sh $1
bash $d/04save.sh $1
