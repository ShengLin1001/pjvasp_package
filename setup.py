from setuptools import setup, find_packages

setup(
name='myalloy-pkg',
version='0.1',
author='B. Yin',
url='#',
description='Package_for_Computational_Metallurgy',
license='#',
packages=find_packages(),
install_requires=['numpy', 'scipy', 'sympy',  'statsmodels', 'ase', 'ovito'],
)