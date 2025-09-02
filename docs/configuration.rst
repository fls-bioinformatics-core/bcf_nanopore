=================================
bcf_nanpore configuration options
=================================

The ``bcf_nanopore`` package can be configured for the local
system by editing the ``bcf_nanopore.ini`` file.

An initial version of this file can be created by copying the
``.sample`` version from the ``config`` subdirectory of the
installed package.

Running ``bcf_nanopre config`` displays the current configuration
and its source.

-----------------------------
General configuration options
-----------------------------

* ``default_runner``: default job runner to use
* ``permissions``: Linux permissions to apply to copied and
  created data (e.g. ``ug+rwX,o-rwX``)
* ``group``: file system group to assign to copied and
  created data

-----------
Job runners
-----------

Define custom runners for certain types of job in this section,
otherwise runners will default to the ``general.default_runner``
specification.

* ``rsync``: runner for large and/or long-running data transfer
  operations (e.g. fetching raw data files)

-------------------
Reporting templates
-------------------

Define custom reporting templates for use by the ``report``
command, for example:

::

   my_report = id,user,#samples,organism,application,PI

Templates are named (``my_report``) and have a list of fields
(separated by commas) to report when the template is used.

Available fields:

- ``name``: project name
- ``id``: project ID
- ``datestamp``: project datestamp
- ``platform``: platform name
- ``user``: associated users
- ``pi``: associated PIs
- ``nsamples``: number of samples
- ``#samples``: alias for ``nsamples``
- ``samples``: comma-separated list of sample names
- ``sample_names``: alias for ``samples``
- ``primary_data``: path to primary data
- ``analysis_dir``: path to the analysis directory
- ``comments``: associated comments
- ``null``: empty value

NB a blank field name is the same as ``null``.
