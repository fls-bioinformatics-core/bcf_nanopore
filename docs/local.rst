=============================================
Local naming and structure conventions at UoM
=============================================

This is an overview of the local conventions adopted at the local
sequencing facility for the naming and structure of PromethION
outputs.

---------------------
Project names and IDs
---------------------

The top-level project folders are named:

::

  PromethION_Project_<NUMBER>_<USER>

For example ``PromethION_Project_002_PerGynt``.

Project identifiers (IDs) take the form:

::

   PROMETHION#<NUMBER>

For example ``PROMETHION#002``

-----------------------------------------------
Output folder structures and naming conventions
-----------------------------------------------

The following naming conventions have been adopted for the
folders under the top-level project directory:

* **Run folders** are named:

  ::

     <SampleRange>_<DDMMYYYY>

  For example, for samples ``PG1`` to ``PG8`` on 26th March 2024
  this would be ``PG1-8_26032024``.

* **Sample pool folders** within a run are named using a loose
  convention which indicates the samples within the pool.

  For example: ``PG1-2`` might be used for a pool containing
  samples ``PG1`` and ``PG2``.

  If sequencing is repeated for a pool then the folders for the
  repeats will have an additional number appended to the name to
  indicate they are repeats (for example ``PG1-2_2``).

(There are currently no conventions for the names or locations of
basecalling output folders outside of flow cell folders.)

-----------------
External metadata
-----------------

While some information can be acquired automatically from folder
names and report files, the following information is expected to
be supplied from the sequencing facility:

* Names of users and PIs
* Application (e.g. "methylation study")
* Organism(s)
* Indexing information relating sample names to flow cells and
  barcodes
