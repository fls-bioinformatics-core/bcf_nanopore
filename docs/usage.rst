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

The new analysis directory will be called ``<PROJECT_NAME>_analysis``,
and will contain a top-level metadata file ``project.info`` along with
subdirectories for each run found in the source data.

Each run subdirectory name will be prefixed with a number (to indicate
the order in which the runs were found), and will contain the following
following files:

* ``flowcell_basecalls.info``: information about each of the flow cell
  and basecalling directories found in the run
* ``sample.info``: a placeholder file for information sample names and
  the associated barcodes and flow cell IDs (this must be filled in
  manually once the analysis directory has been created)
* HTML report files from the basecalling (if found in the source data)

The following command line options are compulsory:

* ``--user USER``: user name
* ``--pi PI``: principal investigator
* ``--application APPLICATION``: application
* ``--organism ORGANISM``: organism

and this information is used to populate the top-level metadata file when
the analysis directory is created.

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

------------
``metadata``
------------

Usage:

::

   bcf_nanopore metdata [--json] REPORT_FILE

Extracts and prints JSON data extracted from ``REPORT_FILE``
(which must either a HTML or JSON report from the PromethION
basecaller.

If the ``--json`` option is supplied then the command dumps
the entire JSON data from the file as-is; otherwise it only
extracts and prints the specific metadata items derived from
the file in question.

(This command is intended to help with manually locating and
identifying useful metadata items within the report data.)
