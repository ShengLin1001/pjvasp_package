#!/bin/bash

head -A 10 POTCAR
cat -n POSCAR

# atomnum+3
grep -A 5 "TOTAL-FORCE" OUTCAR

grep "in kB" OUTCAR