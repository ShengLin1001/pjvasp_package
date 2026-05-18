#!/bin/bash

# VASP
head -A 10 POTCAR
cat -n POSCAR
grep -A 5 "TOTAL-FORCE" OUTCAR    # atomnum+3
grep "in kB" OUTCAR

# VS code
ctrl + shift + P # open settings.json
ctrl + /   # comment or uncomment line

# codex
ctrl + N # new chat

# internet
curl -I https://www.baidu.com
curl ipinfo.io
curl cip.cc