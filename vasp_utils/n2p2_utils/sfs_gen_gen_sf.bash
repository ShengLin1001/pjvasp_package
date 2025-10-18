#!/bin/bash
# Generated SF requests
SFGEN='basic_SF.py'

mypython=/home/tliu/usr/bin/python3.8/bin/python3.8
#/home/lli/software/python37/bin/python3.7

cp input.nn ../active_sf_0/

###   element -c cutoff -n number_of_cutoff(2n+1) -z specific_parameter


$mypython $SFGEN -G G2 -e Mg -c 3.3 -n 5    > sflist
$mypython $SFGEN -G G2 -e Mg -c 4.6 -n 5  >> sflist
$mypython $SFGEN -G G2 -e Mg -c 5.6 -n 5  >> sflist
$mypython $SFGEN -G G2 -e Mg -c 6.2 -n 5    >> sflist
$mypython $SFGEN -G G3 -e Mg -c 3.3 -n 5 -z 1,4,16  >> sflist
#cat sflist >> input.nn
#cp input.nn ../active_sf_0/input.nn.all
