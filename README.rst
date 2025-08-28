bcf_nanopore
============

.. image:: https://github.com/fls-bioinformatics-core/bcf_nanopore/workflows/Python%20CI/badge.svg
   :target: https://github.com/fls-bioinformatics-core/bcf_nanopore/actions?query=workflow%3A%22Python+CI%22

Overview
--------

Provides Python utilities to help with the management of sequencing
data from Oxford Nanopore Technologies (ONT) instruments (specifically
PromethION data) within the Bioinformatics Core Facility (BCF) at the
University of Manchester.

Installation
------------

It is recommended to install the code by first ``git cloning`` the
repository, and then using ``pip`` to install into a Python virtual
environment using ``pip``, e.g.

::

   $ git clone https://github.com/fls-bioinformatics-core/bcf_nanopore.git
   $ virtualenv -p python3 ./venv
   $ source ./venv/bin/activate
   $ pip install -r ./bcf_nanopore/requirements.txt
   $ pip install ./bcf_nanopore/
    

*NB please use these instructions only as a starting point*

Configuration
-------------

The installed package can be configured via an ``.ini``-style file
called ``bcf_nanopore.ini`` placed either in the current directory
or in the ``config`` directory of the package installation.

.. note::

   A sample version of the config file called
   ``bcf_nanopore.ini.sample``can be found in the ``config``
   directory.

Usage
-----

The package provides a single utility called ``bcf_nanopore`` which
in turn provides a set of subcommands.

To get general help run:

::

   $ bcf_nanopore -h

The main commands are:

* ``config``: print current configuration information
* ``info``: reports information about a PromethION project directory
  containing outputs from one or more runs of the sequencer
* ``setup``: creates an "analysis" directory for a PromethION
  project
* ``fetch``: copies a subset of data from the sequencer output for
  use with downstream analysis
* ``report``: reports the metadata from a PromethION project analysis
  directory
