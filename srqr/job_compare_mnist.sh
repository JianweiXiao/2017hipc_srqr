#!/bin/bash
# m: number of rows of input matrix
# n: number of columns of input matrix
# k: target rank
# nb_max: block size
# p: oversampling size
# g: user defined parameter for the extra swaps stage
# d: number of rows of random matrix used in extra swaps stage
# counter: select test matrix
m=60000
n=784
k=300
nb_max=64
p=10
g=5.0D+0
d=10
counter=2

for k in 50 100 150 200 250 300 350 400 450 500 550 600
do	./test_srqr $m $n $k $nb_max $p $g $d $counter
done