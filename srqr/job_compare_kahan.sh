#!/bin/bash
# m: number of rows of input matrix
# n: number of columns of input matrix
# k: target rank
# nb_max: block size
# p: oversampling size
# g: user defined parameter for the extra swaps stage
# d: number of rows of random matrix used in extra swaps stage
# counter: select test matrix
m=96
n=96
k=95
nb_max=64
p=10
g=5.0D+0
d=10
counter=3

./test_srqr $m $n $k $nb_max $p $g $d $counter

m=192
n=192
k=191
nb_max=64
p=10
g=5.0D+0
d=10
counter=3

./test_srqr $m $n $k $nb_max $p $g $d $counter

m=384
n=384
k=383
nb_max=64
p=10
g=5.0D+0
d=10
counter=3

./test_srqr $m $n $k $nb_max $p $g $d $counter