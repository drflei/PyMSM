"""Setup script for PyMSM package."""
import os
from setuptools import setup, find_packages

# Get requirements for installation
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = os.path.join(lib_folder, 'requirements.txt')
install_requires = []
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

# Read long description from README
readme_path = os.path.join(lib_folder, 'README.md')
long_description = ''
if os.path.isfile(readme_path):
    with open(readme_path, encoding='utf-8') as f:
        long_description = f.read()

setup(
    name='pymsm',
    version='0.1.1',
    description='Python library for calculating the geomagnetic rigidity cutoff',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='drflei',
    author_email='',
    url='https://github.com/drflei/PyMSM',
    license='LGPL',
    py_modules=['pymsm'],
    package_data={'': ['MAPS/**/*.AVG']},
    install_requires=install_requires,
    include_package_data=True,
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords='geomagnetic rigidity cutoff magnetosphere space physics',
    project_urls={
        'Source': 'https://github.com/drflei/PyMSM',
        'Bug Reports': 'https://github.com/drflei/PyMSM/issues',
    },
)
