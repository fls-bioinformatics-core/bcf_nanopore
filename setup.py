"""Description

Setup script to install bcf_nanopore

Copyright (C) University of Manchester 2024-2025 Peter Briggs

"""

# Hack to acquire all scripts that we want to
# install into 'bin'
from glob import glob
scripts = [ "bin/bcf_nanopore" ]
for pattern in ('bin/*.py','bin/*.sh',):
    scripts.extend(glob(pattern))

# Installation requirements
install_requires = ['genomics-bcftbx',
                    'auto-process-ngs']

# If we're on ReadTheDocs then we can reduce this
# to a smaller set (to avoid build timeouts)
import os
if os.environ.get("READTHEDOCS") == "True":
    install_requires = []

# Setup for installation etc
from setuptools import setup
import bcf_nanopore
setup(name = "bcf_nanopore",
      version = bcf_nanopore.get_version(),
      description = "Utilities to handle Oxford Nanopore data",
      long_description = """Utilities to help with handling data from
      Oxford Nanopore Technologies sequencing instruments""",
      url = "https://github.com/fls-bioinformatics-core/bcf_nanopore",
      maintainer = "Peter Briggs",
      maintainer_email = "peter.briggs@manchester.ac.uk",
      packages = [ "bcf_nanopore",
                   "bcf_nanopore.nanopore", ],
      license = 'AFL-3',
      # Pull in dependencies
      install_requires = install_requires,
      # Enable 'python setup.py test'
      test_suite="nose.collector",
      tests_require=[ "nose" ],
      # Scripts
      scripts = scripts,
      # Sample configuration file
      data_files = [("config", ["config/bcf_nanopore.ini.sample"]),],
      classifiers=[
          "Development Status :: 2 - Pre-Alpha",
          "Environment :: Console",
          "Intended Audience :: End Users/Desktop",
          "Intended Audience :: Science/Research",
          "Intended Audience :: Developers",
          "License :: OSI Approved :: Academic Free License (AFL)",
          "Operating System :: POSIX :: Linux",
          "Operating System :: MacOS",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Programming Language :: Python :: 3",
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
      ],
      include_package_data=True,
      zip_safe = False)
