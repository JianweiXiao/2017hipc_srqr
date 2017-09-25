m=1000
n=1000
k=1000
nprow=2
npcol=2
nb_alg=64
np_oversampling=10
nb_dist=64

for k in 1000
do
    mpirun -n 4 ./xtimer_qr -m $m -n $n -k $k -p $nprow -q $npcol -nb_alg $nb_alg -np_oversampling $np_oversampling -nb_dist $nb_dist
done
