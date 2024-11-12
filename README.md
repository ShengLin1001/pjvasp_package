# pjvasp_package

A package designed for modeling and postprocessing metal systems, with a special focus on thin film applications.

## Install this package:

```shell
pip install -r requirements.txt
pip install git+https://github.com/BinglunYin/myalloy_package.git
###########
pip install git+https://github.com/ShengLin1001/pjvasp_package.git
# or, download the repository.
pip install -e  .
###########
pip uninstall mymetal-pkg
```

## Structure: The documentation is currently under construction using sphinx.

```shell
./mymetal/
├── build
│   └── film
├── calculate
│   ├── calenergy
│   ├── calmath
│   ├── calmechanics
│   ├── calmismatch
│   ├── calqm
│   └── electronic_structure
├── example
│   ├── test-cut
│   ├── test-generate-bulk
│   ├── test-hetbuilder-fixatom
│   ├── test-hydroxylated
│   ├── test-hydroxylated-custom
│   ├── test-post
│   ├── test-stack
│   ├── test-stretch
│   └── test-surface-energy
├── io
│   ├── post
│   └── vasp.py
├── ml
│   ├── confusionmatrix.py
│   ├── dataset.py
│   ├── model.py
│   └── plot.py
├── post
│   ├── newmain.py
│   └── oldmain.py
└── universial
    ├── atom
    ├── check
    ├── data
    ├── index
    ├── matrix
    ├── plot
    ├── print
    └── search

33 directories, 7 files
```

Each module and function includes a docstring. If you have any questions, please refer to the source code or use the help() function.

I am working on generating the documentation directly from the docstrings. The /mymetal/example/ directory contains all the special functions I have developed.

## Additional requirements for ml, heterostructure.

```shell
############## jupyter
# iprPy

############## ml package
# torch
# torchvision
# scikit-learn

############## find heterostructure
# hetbuilder
```

