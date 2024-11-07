#!/bin/bash
#tree ./mymetal/ -I "__pycache__|*.pyc" > structure.txt
tree -L 2 ./mymetal/ -I "__pycache__|*.pyc|__init__.py" > outline.txt

