#!/bin/bash

mypy=/home/tliu/usr/bin/python3.8/bin/python3.8
#/home/lli/software/python37/bin/python3.7

#prepare for lammps
cp input.nn result/post/test/lammps/nnp-data/input.nn
cp scaling.data result/post/test/lammps/nnp-data/scaling.data

#get post
cp input.data result/post/input.data
cp input.nn result/post/input.nn
cp learning-curve.out  result/post/learning-curve.out

cd result/post
#$mypy test_learning_curve_slow.py 
$mypy train_learning_curve_fast.py 

# learning_curve.py


echo 'success'
