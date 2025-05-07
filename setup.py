import os
import sys
import subprocess
from setuptools import setup

requirements=[
    'attrs>=23.2.0',
    'demes>=0.2.3',
    'msprime>=1.3.1',
    'newick>=1.9.0',
    'numpy>=1.26.4',
    'PyYAML>=6.0.1',
    'referencing>=0.35.1',
    'rpds-py>=0.18.1',
    'ruamel.yaml>=0.18.6',
    'ruamel.yaml.clib>=0.2.8',
    'svgwrite>=1.4.3',
    'tskit>=0.5.6'
]

setup(
    name='adnator',
    install_requires=requirements,
    version='1.3.2',
    author='Jazeps Medina Tretmanis',
    author_email='jaz.medtre@gmail.com',
    description='A realistic-ish aDNA simulator',
    url='https://github.com/Jazpy/adnator',
    classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: POSIX :: Linux"],
    packages=['adnator'],
    package_dir={'adnator': 'src/adnator'},
    package_data={'adnator': ['*.cpp', '*makefile']},
    include_package_data=True,
)
