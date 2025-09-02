========================================
Background: PromethION Sequencer Outputs
========================================

------------
Introduction
------------

This document aims to give an overview of the structure of data
produced from runs of the Oxford Nanopore Technologies (ONT)
PromethION instrument.

--------
Overview
--------

Each run of the PromethION instrument involves sequencing a pool
of barcoded samples which have been loaded into one or more flow
cells. Individual samples can be identified using the combination
of a flow cell ID plus a barcode ID. The raw data from the run is
encoded into POD5 format.

Optionally basecalling can be run in real time while the data are
being collected, generating FASTQ and/or BAM files. It is also
possible to (re)run the basecalling on previously collected POD5
data.

A project is made up of data from multiple runs.

-----------------
Glossary of terms
-----------------
  
Within this document the following terminology is used to refer
to various conceptual elements of the data structure:

* PROJECT: a set of runs relating to the same experiment, with
  samples run for one or more applications
* RUN: a run of the PromethION instrument on a pool of samples
* POOL: a set of samples run on the same flow cell within a run
* FLOWCELL: a physical flow cell used in a run, identified by
  a unique ID (e.g. "PAW15685")
* BARCODE: a barcode sequence attached to a sample, identified
  by an integer number (1-24)
* SAMPLE: a named physical sample provided by a user (e.g.
  "PG1")
* APPLICATION: the type of experiment (e.g. "methylation study")
* BASECALLS DIRECTORY: a directory or folder containing base
  calls data

This terminology differs slightly from that seen elsewhere in
PromethION documentation:

* A run in our terminology may be referred to as an "experiment"
* A pool in our terminology may be referred to as a "sample"

-----------------------------------
PromethION sequencer output folders
-----------------------------------

************************
Projects, runs and pools
************************

Within the local sequencing facility, a single top-level data
folder will conventionally hold the outputs from one or more run
associated with the same project.

The PromethION imposes very few restrictions on the naming or
structure of folders within each project data folder, however
within the sequencing facility this top-level data folder should
have the following general structure:

::

   <PROJECT>/[.../][<RUN>/]<POOL>/<FLOWCELL>/bam_pass/<BARCODE>/...

For example:

::

   PromethION_Project_002_PerGynt/PG21-30_12062024/PG21-22/20240612_0123_1A_PAW12345_678ab90c/bam_pass/barcode01/...

where:

* PROJECT = "PromethION_Project_002_PerGynt"
* RUN = "PG21-30_12062024"
* POOL = "PG21-22"
* FLOWCELL = "20240612_0123_1A_PAW12345_678ab90c"
* BARCODE = "barcode01"

Variations on this structure include:

* Missing RUN folder level may be absent (e.g. ``PromethION_Project_002_PerGynt/PG21-22/20240612_0123_1A_PAW12345_678ab90c/bam_pass/barcode01/...``);
* Additional arbitrary folder levels present between the PROJECT
  and RUN levels.
  
Within each run folder (or under the top-level project folder,
if no run folders are present) there may then be one or more
further subfolders which correspond to pools of samples which
have been run on the same flow cell. These pool subfolders will
have names like e.g. "PG1-2".

While the concepts of runs and pools reflect the physical processes
of running the experiments, they are may not necessarily be
reflected in the final data folder structure.

Ultimately the two "core units of outputs" are:

* Flow cell folders
* Basecalling folders

The structure and contents of these two types of output folder are
explained in more detail in the following sections.

************************
Flow cell output folders
************************

Each flow cell folder holds outputs from a run of the PromethION
instrument for a pool of samples run on a single physical flow
cell.

Flow cell folders have a machine-generated name of the form:

::

   <YYYYMMDD>_<HHMM>_<??>_<FlowCellID>_<???>

(e.g. "20240328_0927_1D_PAW15685_dc41c950").

Each flow cell folder has a standard set of subfolders containing
the data produced by the instrument and the basecalling software:

* ``pod5_pass`` and ``pod5_fail`` folders with the raw POD5 data
  should always be present
* Other ``*_pass`` and ``*_fail`` folders hold the outputs from the
  basecalling software in different formats (e.g. ``bam_pass`` for
  BAM files)

The ``*_pass`` folders are in turn are subdivided into folders for
each barcode number (e.g. ``barcode04``).

As each physical sample is run in a specific flow cell with a
specific barcode, it is possible to use these to locate the data for
each sample. For example: if sample ``PG15`` was run in the flow
cell with ID "PAW16585" and associated barcode number 17, then the
corresponding BAM files should be under:

::

   PromethION_Project_002_PerGynt/PG9-20_22042024/PG15-16/20240422_1039_1D_PAW16585_eb5adc49/bam_pass/barcode17

Additionally a number of report files output from basecalling may
be present in the flow cell folder; these are described in a
separate section below.

**************************
Basecalling output folders
**************************

It is possible to (re)run the basecalling software on the POD5
data from the flow cell folders.

The location and names of the basecalling output folders are
completely arbitrary, and the structure of the data in these folders
differs from that of the flow cell output folders: it consists of
two subfolders ``pass`` and ``fail``, each of which contains a
set of ``barcode*`` subfolders (for each barcode).

For example:

::

   PromethION_Project_002_PerGynt/PG21-30_13052024/Rebasecalling/PG21-22/pass/barcode01

These barcode subfolders will contain both BAM and FASTQ files
output by the basecaller.

************************
Basecalling report files
************************

A number of report files produced by the basecalling software may
also be present in the flow cell and basecalling output folders.

These will be named:

::

   report_<FlowCellID>.[html|json|md]

These files are not guaranteed to be present (for example if the
real time basecalling was interrupted, or if the appropriate
options were not selected when setting up the basecalling).

The files appear to contain overlapping sets of metadata and
other information, which can be extracted with some effort to
provide information about various aspects of the basecalling.

The HTML report can be viewed in an appropriate web browser.
