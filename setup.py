# -*- coding: utf-8 -*-
"""
Created on May 2024

    Author: P. Ferraiuolo
"""
from setuptools import setup, find_packages

setup(
    name='ggcas',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'astropy',
        'astroquery',
        'scipy'
    ],
    author='Pietro Ferraiuolo',
    author_email='pietro.ferraiuolo@inaf.it',
    description='Software to gather and analyze astrometry and photometry data from Gaia',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/my_package',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
