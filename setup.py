import sys
from pathlib import Path

from setuptools import setup, find_packages


def load_requirements(path):
    requirements = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        package = line.split("#", 1)[0].strip()
        if not package:
            continue
        if package.lower() == "pyside2" and sys.version_info >= (3, 11):
            continue
        if package not in requirements:
            requirements.append(package)
    return requirements


required_packages = load_requirements("requirements.txt")

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
