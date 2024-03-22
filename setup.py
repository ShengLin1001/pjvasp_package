from setuptools import setup, find_packages

setup(
    name='mymetal-pkg',
    version='0.1',
    author='J. Pei',
    author_email='zju_pj@163.com',
    url='https://github.com/ShengLin1001/pjvasp_package.git',
    description='Package_for_Computational_Metallurgy',
    license='#',
    packages=find_packages(),
    install_requires=['numpy', 
                      'scipy', 
                      'sympy',  
                      'statsmodels', 
                      'ase == 3.22.1', 
                      'ovito'],
)
