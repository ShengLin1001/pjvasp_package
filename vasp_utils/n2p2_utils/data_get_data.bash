#!/bin/bash

mypython=/home/tliu/usr/bin/python3.8/bin/python3.8
#/home/lli/software/python37/bin/python3.7

$mypython 'read.py'
$mypython 'input.py'
cp input.data  ../train/
cp input.data  ../active_sf_0/
cp forces.npy ../active_sf_0/ 
cp energy.npy ../active_sf_0/ 
cd vasp_test
cp -r * ../../train/result/post/test/vasp/
