from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required_packages = f.read().splitlines()

setup(
    name='mymetal-pkg',
    version='1.0.0',
    author='J. Pei',
    author_email='zju_pj@163.com',
    url='https://github.com/ShengLin1001/pjvasp_package.git',
    description='A comprehensive package for computational metallurgy, providing tools for material property calculations, structure modeling, and analysis.',
    license='MIT',
    packages=find_packages(),
    install_requires=required_packages
)
