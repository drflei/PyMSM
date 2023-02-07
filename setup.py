from distutils.core import setup
from setuptools import find_packages
import os

# get requirements for installation
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'
install_requires = [] # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

setup(
    name='pymsm',
    packages=find_packages(','),
    
    # package_data={"":["DLRISOmodel/neutronMonitorData/*.dat","DLRISOmodel/neutronMonitorData/*.pkl"]},
    include_package_data=True,
    version='0.1.1',
    description='Python library for calculating the geomagnetic rigidity cutoff',
    author='drflei',
    license='LGPL',
    install_requires=install_requires,
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='test',
)
