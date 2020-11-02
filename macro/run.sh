export MAT='G4_BRAIN_ICRP'

for MAT in 'G4_BRAIN_ICRP' 'Cu63' 'Sc45' 'Y89' 'F19'
do
envsubst '$MAT' < run.mac > run_$MAT.mac
envsubst '$MAT' < /mnt/data1/cartechini/runGeantMAT_1e8.sh > /mnt/data1/cartechini/runGeant${MAT}_1e8.sh
qsub /mnt/data1/cartechini/runGeant${MAT}_1e8.sh -q move-it
done
