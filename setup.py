from setuptools import setup
from codecs import open
import os
import re

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),'sereader','__version__.py')

with open(version_file, "r") as f:
    lines = f.read()
    sereader_version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", lines, re.M).group(1)

setup(
    name='sereader',
    version=sereader_version,
    author='Kristof Rozgonyi',
    author_email='rstofi@gmail.com',
    description='Python package to read serial binary image or spectra (.ser) created with FEI, Tecnai or Titan microscopes',
    install_requires=requirements
    )