#!/bin/bash

mypy=/home/tliu/usr/bin/python3.8/bin/python3.8
#/home/lli/software/python37/bin/python3.7

#echo which epoch you need
#read epoch_number
##read from script
#if [ ! $epoch_number]; then
#epoch_number=$1
#fi

read epoch_number < result/post/epoch.txt   

#get nnp-data
cp weights.012.${epoch_number}.out result/nnp-data/weights.012.data
#cp weights.079.${epoch_number}.out result/nnp-data/weights.079.data
cp input.nn result/nnp-data/input.nn
cp scaling.data result/nnp-data/scaling.data

#get post
cp trainpoints.${epoch_number}.out result/post/trainpoints.data
cp trainforces.${epoch_number}.out result/post/trainforces.data

#prepare nnp to eval test*
cp weights.012.${epoch_number}.out result/post/test/lammps/nnp-data/weights.012.data
#cp weights.079.${epoch_number}.out result/post/test/lammps/nnp-data/weights.079.data

#get test* of this epoch
cd result/post/test
./eval_test.bash
cd ../
printf $epoch_number > epoch.txt
$mypy compare.py

echo 'success'

