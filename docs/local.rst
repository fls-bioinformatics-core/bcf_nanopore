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

Under the top-level project folder, it is expected that there
will be one or more run folders, corresponding to batches of
samples run at different times.

There is no official naming convention used for the run folders
and their internal structure is also not formally defined.
Where a convention is used, it is likely that the run name
will either replicate the project name (for single-run projects)
or else use some indication of the samples included in the run,
for example:

::

    <SampleRange>_<DDMMYYYY>

For example, for samples ``PG1`` to ``PG8`` on 26th March 2024
this would be ``PG1-8_26032024``.

It is possible that additional intermediate folders may be
present between the run folder and the flow cell folders, for
example to indicate pools of samples run together. In this case
the "pool" folders use a loose naming convention which indicates
the samples within the pool.

For example: ``PG1-2`` might be used for a pool containing
samples ``PG1`` and ``PG2``. Within this convention, if
sequencing is repeated for a pool then the folders for the
repeats will have an additional number appended to the name to
indicate they are repeats (for example ``PG1-2_2``).

However the model of the data used in the Python package makes
no assumptions about the presence or naming of these intermediate
sample pool folders (or any other intermediate folders).

Similarly, the names and locations of any basecalling output
folders are assumed to be arbitrary.

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
