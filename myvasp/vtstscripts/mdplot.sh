#!/bin/sh
if [ "$1" == "-h" ]; then
    echo "usage: mdplot.sh [OUTCAR]"
    echo "       outputs a gnuplot script that plots the energy and"
    echo "       temperature vs time"
    echo 
    echo "example: mdplot.sh | gnuplot -p"
    echo
    exit 0
fi
if [ -z "$1" ]
then
    fn=OUTCAR
else
    fn=$1
fi

timestep=`awk '/POTIM/ && $3 ~ /[0-9].*/ {print $3}' $fn`

title=$(basename $(pwd))
echo "set title \"$title\""
echo 'set ytics nomirror'
echo 'set xlabel "Time (fs)"'
echo 'set y2tics auto'
echo 'set ylabel "Temperature (K)"'
echo 'set y2label "Energy (eV)"'
echo 'plot "-" w l t "Temperature (K)", "-" w l axis x1y2 t "Total Energy"'
awk "/EKIN/ {n+=1*$timestep;print n, \$7}" $fn
echo 'e'
awk "/total energy   ETOTAL =/ {n+=1*$timestep;print n, \$5}" $fn
#echo 'e'
#awk "/EKIN/ {n+=1;sum += \$7} END {print sum/n}" $fn
