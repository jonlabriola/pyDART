from setuptools import setup, find_packages, Extension

import os

import numpy as np

#import pycaps._version as vers

np_loc = os.path.dirname(np.__file__)


setup(
    name='pydart',
    version='1.0',
    description='DART processing code and visualization for the CM1 model.',
    packages=find_packages(),
    package_data={"pydart": ['data/DART_obs.csv']},
    include_package_data=True,

    author="Jonathan Labriola",
    author_email="jonathan.labriola@noaa.gov",
)
