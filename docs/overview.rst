===================================================================
bcf_nanpore: utilities for managing Oxford Nanopore sequencing runs
===================================================================

--------
Overview
--------

``bcf_nanopore`` provides a set of Python utilities for managing the
outputs from Oxford Nanopore Technology (ONT) sequencing runs within
the Bioinformatics Core Facility (BCF) at the University of
Manchester.

NB currently only outputs from PromethION sequencers are supported.

The utilities support the following workflow:

1. ONT sequencing instrument outputs data from multiple runs to a
   **project** directory on a shared drive
2. An **analysis directory** is created for the project, to store
   metadata and results from downstream analyses
3. A subset of data files from the project (typically BAMs) are
   copied to a separate **data directory** for use in downstream
   analyses

There also are functions to support examining the contents of (and
metadata associated with) project directories, and reporting information
stored in analysis directories (e.g. for book-keeping and auditing).

It is also possible to configure the package for the local systems.

--------
Contents
--------

* `Background information on the structure of PromethION sequencer outputs <background.rst>`_
* `Configuration <configuration.rst>`_
