#!/bin/bash
cp ../../../y_full_relax/* ./
bash clean_up.sh
bash sed_incar.sh ISTART 0
bash sed_incar.sh NSW 100
bash sed_incar.sh ISYM 1
