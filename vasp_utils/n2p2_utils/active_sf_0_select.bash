#!/bin/bash

mypython=/home/tliu/usr/bin/python3.8/bin/python3.8
#/home/lli/software/python37/bin/python3.7

#select number of SFs
$mypython collect_sf.py
$mypython select_feat.py 48

cp input.nn  ../train/
