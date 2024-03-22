# pjvasp_package

Package for DFT calculations

# Install this package:

## --user for all virtual environment

```shell
pip3 install --upgrade --user git+https://github.com/ShengLin1001/pjvasp_package.git
```

## Recommended
Firstly, you need to download .zip and unzip it (Insteading, just cloning using desktop GitHub!)

Secondly, enter the family directory, and run the following bash command

```shell
pip3 install --user -e  .  
```

# some fixed problems

## when i run the test-post/*.py, i found the output file has been generated in the root directory of pjvasp_package

we should add fllowing python commands

```shell
import os
# get the path of those script
script_dir = os.path.dirname(__file__)
# set the working path
os.chdir(script_dir)
```

# some fixing problems


