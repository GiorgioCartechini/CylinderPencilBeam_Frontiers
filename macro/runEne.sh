export ENE=0

for ENE in 5 20
do
envsubst '$ENE' < runEne.mac > run_$ENE.mac
envsubst '$ENE' < /mnt/data1/cartechini/runGeantENE_2e9.sh > /mnt/data1/cartechini/runGeant${ENE}_2e9.sh
qsub /mnt/data1/cartechini/runGeant${ENE}_2e9.sh -q move-it
done
