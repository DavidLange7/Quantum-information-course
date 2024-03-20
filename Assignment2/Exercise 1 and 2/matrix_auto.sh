#!/bin/bash
gfortran matrix.f90 -o matrix
touch one.txt
rm one.txt
touch one.txt
touch two.txt
rm two.txt
touch two.txt
touch three.txt
rm three.txt
touch three.txt
./matrix <<!
0
1
1
1
1
!

for n_temp in {15..50}; do
./matrix <<!
0
$n_temp
$n_temp
$n_temp
$n_temp
!
done

for n_glob in {2..30}; do
let temp="50*n_glob"
./matrix <<!
0
$temp
$temp
$temp
$temp
!
done
