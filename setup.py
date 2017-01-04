from setuptools import setup, find_packages
import sys, os

VERSION = '1.2.0'

setup(name='garmire_SSrGE',
      version=VERSION,
      description="compute SNV from RNA-seq following GATK recommendations",
      long_description="""""",
      classifiers=[],
      keywords='',
      author='Olivier Poirion (PhD)',
      author_email='opoirion@hawaii.edu',
      url='',
      license='MIT',
      packages=find_packages(exclude=['examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=['numpy',
                        'scipy',
                        'scikit-learn==0.18',
                        'tabulate'],
      )
