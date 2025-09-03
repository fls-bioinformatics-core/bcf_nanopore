=================
bcf_nanpore usage
=================

-------
General
-------

There is a single command ``bcf_nanopore`` which implements various
functions as subcommands.

The general usage format is:

::

   bcf_nanopore CMD [options]

For example:

::

   bcf_nanopore fetch /mnt/data/PromethION_Project_004

The documentation below gives the basic information about each of
the commands; use ``bcf_nanopore CMD -h`` to see all options.

----------
``config``
----------

Usage:

::

   bcf_nanopore config

Display the current configuration information.

--------
``info``
--------

Usage:

::

   bcf_nanopore info <PROJECT>

Analyses a PromethION project data folder and prints information
about flow cells, basecalling directories and extracted metadata.

---------
``setup``
---------

Usage:

::

   bcf_nanopore setup <PROJECT> <DEST> ...

Creates an analysis directory under ``DEST`` for the PromethION
project ``PROJECT``.

The new analysis directory will be called ``<PROJECT_NAME>_analysis``.

The following command line options are compulsory:

* ``--user USER``: user name
* ``--pi PI``: principal investigator
* ``--application APPLICATION``: application
* ``--organism ORGANISM``: organism

This information is added to the metadata file when the analysis
directory is created.

Optionally a CSV file which list the samples along with their
associated flow cells and barcodes can also be supplied via the
``--samples_csv`` option. The file can have two fields (sample
name and barcode) or three (sample name, barcode and flow cell ID).

For example:

::

   #Sample,Barcode,Flowcell
   SMPL_A1,BP03,PBC31213
   SMPL_A2,BP04,PBC31213
   SMPL_B1,BP01,PBC31213
   SMPL_B2,BP02,PBC31213
   
If supplied then this information is also added to the analysis
directory metadata.

---------
``fetch``
---------

Usage:

::

   bcf_nanopore fetch PROJECT DEST [--files FILETYPES]
   
Copies a subset of the Promethion output data from ``PROJECT``
under ``DEST``.

``PROJECT`` can be a local or remote directory.

By default the subset consists of BAM files from the basecalling
plus any report files found in the source; the ``--files`` option
can be used to specify the set of file types, for example:

::

   --files pod5,fastq,bam

will copy POD5, FASTQ and BAM files, whereas

::

   --files fastq,bam

will only copy FASTQ and BAM files.

----------
``report``
----------

Usage:

::

   bcf_nanopore report ANALYSIS_DIR

Report metadata from ``ANALYSIS_DIR``.
