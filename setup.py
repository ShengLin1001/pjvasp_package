from setuptools import setup, find_packages

# 读取requirements.txt文件中的内容到列表中
with open('requirements.txt') as f:
    required_packages = f.read().splitlines()

setup(
    name='mymetal-pkg',
    version='0.2.0',
    author='J. Pei',
    author_email='zju_pj@163.com',
    url='https://github.com/ShengLin1001/pjvasp_package.git',
    description='Package_for_Computational_Metallurgy',
    license='#',
    packages=find_packages(),
    install_requires=required_packages,  # 使用从requirements.txt读取的依赖项
)
