#!/usr/bin/env python3
#
#     analysis: handle analyses of ONT PromethION data
#     Copyright (C) University of Manchester 2025 Peter Briggs
#

"""
Library with classes for analysing Oxford Nanopore Technologies
PromethION data for BCF.

Provides the following classes:

* ProjectAnalysisDir: analysis directory associated with a project
* ProjectInfo: metadata about a project
* SamplesInfo: metadata about samples in a project
* FlowcellBasecallsInfo: metadata about flow cell basecalls
"""

import os
import re
import shutil
import logging
from pathlib import Path
from auto_process_ngs.metadata import MetadataDict
from auto_process_ngs.utils import get_numbered_subdir
from bcftbx.TabFile import TabFile
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from .nanopore.promethion import ProjectDir
from .nanopore.promethion import get_flow_cell_datestamp
from .utils import MetadataTabFile
from .utils import convert_field_name
from .utils import fmt_value
from .utils import fmt_yes_no

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RE_PROJECT_DIR_NAME = re.compile("^PromethION_Project_([0-9]+)_(.+)$")
RE_RUN_DIR_NAME = re.compile("^([0-9]{3})_(.+)$")

class ProjectAnalysisDir:
    """
    Create and manage directory for analysing PromethION project

    The analysis directory is a place to store metadata about
    a project, and to put any data produced by downstream
    analysis.

    The basic structure looks like:

    PROJECT_analysis
      |
      +-- README
      |   project.info
      |
      +-- 01-Run_name
      |     |
      |     +-- flowcell_basecalls.tsv
      |         report...
      |
      +-- 02-Run_name
      |     |
      |     ...
      |
      +-- ScriptCode
      |
      +-- logs
    
    An analysis directory can be created for a PromethION
    project using the ``create`` method.

    Arguments:
      path (str): path to the analysis directory
    """

    def __init__(self, path):
        # Top-level metadata
        self.path = os.path.abspath(path)
        print(self.path)
        # Load top-level metadata
        self.project_info_file = os.path.join(self.path,
                                              "project.info")
        self.info = ProjectInfo()
        if os.path.exists(self.project_info_file):
            self.info.load(self.project_info_file)
        else:
            logger.warning("%s: no 'project.info' file found" %
                           self.path)
        # Dictionary of run_dirs and their directories
        self.run_dirs = {}
        if self.exists():
            for d in os.listdir(self.path):
                if os.path.isdir(os.path.join(self.path, d)) and RE_RUN_DIR_NAME.match(d):
                    run_name = RE_RUN_DIR_NAME.match(d).group(2)
                    self.run_dirs[run_name] = d

    def create(self, project_dir, user, PI, application, organism):
        """
        Creates a new PromethION project analysis directory

        Requires basic metadata data items (user, PI, application,
        organism).

        Optionally information about samples in the project
        can also be supplied.

        Arguments:
          project_dir (str): path to PromethION project
            directory
          user (str): name of user
          PI (str): name of PI (principal investigator)
          application (str): application type
          organism (str): associated organism(s)
        """
        if self.exists():
            raise OSError("%s: already exists" % self.path)
        os.makedirs(self.path)
        # Store top-level metadata
        metadata = {
            'user': user,
            'PI': PI,
            'application': application,
            'organism': organism
        }
        for item in metadata:
            if metadata[item] is not None:
                self.info[item] = str(metadata[item]).strip()
        self.info['data_dir'] = os.path.abspath(project_dir)
        self.info['name'] = os.path.basename(self.info['data_dir'])
        self.info['id'] = self._make_project_id(self.info['name'])
        self.info['platform'] = "promethion"
        self.info.save(filen=self.project_info_file)
        # Load source project directory
        project = ProjectDir(project_dir)
        # Handle run_dirs
        for run in project.runs:
            # Create a subdirectory for the run
            print(f"-- run '{run.name}'")
            run_dir = get_numbered_subdir(run.name, self.path, full_path=True)
            print(f"   creating directory: '{run_dir}'")
            os.makedirs(run_dir)
            print(f"   populating run directory...")
            # Create a TSV file with flow cell and base calls info
            flow_cell_basecalls_file = os.path.join(run_dir,
                                                    "flowcell_basecalls.tsv")
            fc_file = FlowcellBasecallsInfo()
            # Add information for flow cells in the run
            for fc in run.flow_cells:
                fc_file.add_base_calls(
                    run=run.name,
                    sub_dir=os.path.relpath(fc.path, project.path),
                    flow_cell_id=fc.id,
                    reports=(",".join(fc.report_types)
                             if fc.report_types else "none"),
                    kit=fmt_value(fc.metadata.kit),
                    modifications=("none"
                                   if fc.metadata.modified_basecalling == "Off"
                                   else fmt_value(fc.metadata.modifications)),
                    trim_barcodes=fmt_value(fc.metadata.trim_barcodes),
                    minknow_version=fc.metadata.software_versions["minknow"],
                    basecalling_model=fmt_value(fc.metadata.basecalling_model),
                    file_types=(",".join(fc.file_types)
                                if fc.file_types else "none"))
            # Add information for basecalls dirs in the run
            for bc in run.basecalls_dirs:
                fc_file.add_base_calls(
                    run=run.name,
                    sub_dir=os.path.relpath(bc.path, project.path),
                    flow_cell_id=fmt_value(bc.metadata.flow_cell_id),
                    reports=(",".join(bc.report_types)
                             if bc.report_types else "none"),
                    kit=fmt_value(bc.metadata.kit),
                    modifications=("none"
                                   if bc.metadata.modified_basecalling == "Off"
                                   else fmt_value(bc.metadata.modifications)),
                    trim_barcodes=fmt_value(bc.metadata.trim_barcodes),
                    minknow_version=bc.metadata.software_versions["minknow"],
                    basecalling_model=fmt_value(bc.metadata.basecalling_model),
                    file_types=(",".join(bc.file_types)
                                if bc.file_types else "none"))
            # Save flow cell/base calls info file
            fc_file.save(flow_cell_basecalls_file)
            # Copy in flow cell base calling reports
            for fc in run.flow_cells:
                report = fc.html_report
                if not report:
                    continue
                target = os.path.join(run_dir,
                                      "%s_%s_%s" % (fc.run,
                                                    fc.id,
                                                    os.path.basename(report)))
                shutil.copy(report, target)
            # Copy in other base calling reports
            for bc in project.basecalls_dirs:
                report = bc.html_report
                if not report:
                    continue
                target = os.path.join(run_dir,
                                      "%s_%s_%s_%s" % (bc.parent,
                                                       bc.run,
                                                       bc.metadata.flow_cell_id,
                                                       os.path.basename(report)))
                shutil.copy(report, target)
            # Create a template samples file
            samples_file = os.path.join(run_dir, "samples.tsv")
            with open(samples_file, "wt") as fp:
                fp.write("#Sample\tBarcode\tFlowcell\n")
            # Create a README file for the run
            read_me_file = os.path.join(run_dir, "README")
            with open(read_me_file, "wt") as read_me:
                read_me.write(
                    f"""This is the analysis directory for run '{run.name}' of project '{self.info.name}'.

The following files have been automatically generated:

- 'flowcell_basecalls.tsv': TSV file listing information about flow cells and
  basecalls directories
""")
                if os.path.exists(os.path.join(run_dir, "samples.tsv")):
                    read_me.write(
                        "- 'samples.tsv': TSV file matching sample names to flow cell "
                        "and barcode IDs\n")
                read_me.write(
                    f"- copies of HTML reports from the basecalling runs found in "
                    "the primary data directory (renamed to identify the associated "
                    "locations).")
            # Add run to list of run_dirs
            self.run_dirs[run.name] = os.path.basename(run_dir)
        # Update the 'runs' field in the project info
        self.info['runs'] = ",".join(self.runs)
        self.info.save(filen=self.project_info_file)
        # Make subdirectories
        for subdir in ("logs", "ScriptCode"):
            Path(self.path).joinpath(subdir).mkdir()
        # Create a README file
        read_me_file = os.path.join(self.path, "README")
        with open(read_me_file, 'wt') as read_me:
            read_me.write(
                f"""This is the analysis directory for {os.path.basename(self.info.data_dir)}

The following files and directories have been automatically generated:

- '{os.path.basename(self.project_info_file)}': top-level information about the project
- Subdirectories for each run in the project, named with a leading
  number to indicate order of processing (e.g. '001_{project.runs[0].name}'):
- 'logs': directory for log files;
- 'ScriptCode': directory for custom scripts and code.
""")
        # Report information
        print(self.report_project_summary())

    @property
    def runs(self):
        """
        Return a list of run names for the project

        Returns:
            list: list of run names
        """
        return sorted(self.run_dirs.keys(), key=lambda x: self.run_dirs[x])

    def exists(self):
        """
        Check whether project analysis directory exists
        """
        return os.path.exists(self.path)

    def report(self, mode, fields=None):
        """
        Report information on the project analysis directory

        Reporting mode can be either "project" (summarises
        runs) or "runs" (explicit output for each run in the
        project).

        Format can be either "summary" (a summary of the project)
        or "runs" (tab-separated values on one line in "project"
        mode, or one line per run in "runs" mode).

        Available fields are those that can be fetched using
        the 'get_value' method.

        Note that some fields may not be available
        on the reporting mode (e.g. 'run' is not available
        in 'project' mode).

        A blank field name is the same as 'null'.

        Arguments:
            mode (str): reporting mode ("summary" or "runs")
            fields (str): optional, comma separated list of field
              names to report in "runs" mode

        Returns:
          String: the report text.
        """
        if mode == "summary":
            return self.report_project_summary()
        elif mode == "runs":
            if fields is None:
                raise Exception(f"'fields' must be specified for 'runs' mode")
            return self.report_project_runs(fields=fields)
        else:
            raise Exception(f"{mode}: not a valid report mode")

    def report_project_summary(self):
        """
        Reports a summary report of the project

        The report consists of a title line, a block of
        key-value pairs of project-level metadata, and
        a block reporting run-specific metadata (one run
        per line) for each run.
        """
        field_names = {
            "name": "Project name",
            "id": "Project ID",
            "user": "User",
            "pi": "PI",
            "application": "Application",
            "organism": "Organism",
            "analysis_dir": "Analysis dir"
        }
        output = []
        # Title
        output = [self.info.name, "="*len(self.info.name), ""]
        # Project-level metadata
        for field in ["name", "id", "user", "pi", "application", "organism", "analysis_dir"]:
            output.append("%-16s: %s" % (field_names[field], self.get_value(field)))
        # Runs
        if self.runs:
            output.extend(["",
                           "This project has %s run%s:" % (len(self.runs), "s" if len(self.runs) != 1 else ""),
                           ""])
            for run in self.runs:
                samples = self.get_value("samples", run)
                if samples:
                    # Convert to a list
                    samples = samples.split(",")
                else:
                    samples = []
                nsamples = len(samples)
                line = ["%s:" % run,
                        "%s sample%s%s" % (nsamples if nsamples else "no",
                                           "s" if nsamples != 1 else "",
                                           " (%s)" % ", ".join(samples) if samples else "")]
                output.append("- %s" % "\t".join([str(x) for x in line]))
        else:
            output.extend(["",
                           "No runs detected for this project?"])
        return "\n".join(output)

    def report_project_runs(self, fields):
        """
        Reports runs in a project

        Reports each run in a project as a tab-separated
        set of values (one run per line).
        """
        fields = [f.strip().lower() for f in fields.split(',')]
        output = []
        for run in sorted(self.runs):
            line = []
            for f in fields:
                line.append(self.get_value(f, run))
            output.append("\t".join([str(x) for x in line]))
        return "\n".join(output)

    def get_value(self, field, run=None):
        """
        Get value of a project information field

        Available fields are:

        - name (project name)
        - id (project ID)
        - datestamp (earliest associated flow cell datestamp)
        - nruns (number of runs)
        - #runs (alias for 'nruns')
        - runs (comma-separated list of run names)
        - platform (platform name)
        - user (associated users)
        - pi (associated PIs)
        - run (name of run)
        - run_datestamp (datestamp of run)
        - nsamples (number of samples)
        - #samples (alias for 'nsamples')
        - samples (comma-separated list of sample names)
        - sample_names (alias for 'samples')
        - primary_data (path to primary data)
        - analysis_dir (path to the analysis directory)
        - comments (associated comments)
        - null (empty value)

        Composite fields can be specified using the syntax

        [DELIMITER]:FIELD1+FIELD2...

        where the DELIMITER is used to separate values from
        the specified FIELDs.

        For example:

        [_]:name+run

        will produce values of the form
        'PromethION_Project_001_PerGynt_PG1-2_20240513'

        Arguments:
          field (str): name of the field to retrieve
          run (str): optional, specifies run to get value
            associated with field (for run-specific fields)

        Returns:
          str: value associated with the field
        """
        delimiter = " "
        field = field.strip().lower()
        if field.startswith("["):
            # Handle custom delimiter for composite field
            field_ = field.split(":")
            if len(field_) > 1 and field_[-2].endswith("]"):
                delimiter = ':'.join(field_[:-1])[1:-1]
            else:
                raise ValueError("Bad delimiter specification")
            field = field_[-1]
        if "+" in field:
            # Composite field
            # Recover value for individual fields
            value = []
            for field in field.split("+"):
                value.append(self.get_value(field, run))
            # Return composite value
            return delimiter.join([str(x) for x in value])
        # Single field specified
        fmt_func = fmt_value
        if field == "" or field == "null":
            value = ''
        elif field == "name":
            value = self.info.name
        elif field == "id":
            value = self.info.id
        elif field == "datestamp":
            value = self.datestamp()
        elif field == "run_datestamp":
            if run is None:
                raise KeyError(f"'{field}' field requires a run name")
            value = self.datestamp(run)
        elif field == "run":
            if run is None:
                raise KeyError(f"'{field}' field requires a run name")
            value = run
        elif field == "runs":
            value = ",".join(self.runs)
        elif field == "nruns" or field == "#runs":
            value = len(self.runs)
        elif field == "platform":
            value = self.info.platform
        elif field == "user":
            value = self.info.user
        elif field == "pi":
            value = self.info.PI
        elif field == "application":
            value = self.info.application
        elif field == "organism":
            value = self.info.organism
        elif field == "nsamples" or field == "#samples":
            if run is None:
                runs = self.runs
            elif run in self.runs:
                runs = [run]
            else:
                raise KeyError("%s: unrecognised run name" % run)
            # Get total number of samples across all runs
            value = 0
            for run in runs:
                samples_file = os.path.join(self.path, self.run_dirs[run], "samples.tsv")
                if os.path.exists(samples_file):
                    value += len(SamplesInfo(samples_file))
            def fmt_func(s): return '?' if s == "" else s
        elif field == 'sample_names' or field == 'samples':
            if run is None:
                runs = self.runs
            elif run in self.runs:
                runs = [run]
            else:
                raise KeyError("%s: unrecognised run name" % run)
            sample_names = []
            for run in runs:
                sample_names.extend([s["Sample"] for s in SamplesInfo(os.path.join(self.path,
                                                                                   self.run_dirs[run],
                                                                                  "samples.tsv"))])
            value = ",".join(sample_names)
        elif field == "primary_data":
            value = self.info.data_dir
        elif field == "analysis_dir":
            value = self.path
        elif field == "comments":
            value = self.info.comments
            def fmt_func(s): return '' if s is None else str(s)
        else:
            raise KeyError("%s: unrecognised field" % field)
        return fmt_func(value)

    def datestamp(self, run=None):
        """
        Fetch a datestamp for project or run

        If no datestamp can be located, return None (datestamps
        are derived from flow cell output directories, with the
        earliest datestamp for the project or run being returned).

        Arguments:
            run (str): optional, specifies run to get datestamp for

        Returns:
            str: earliest associated datestamp extracted from
            flow cell directories within project or run
        """
        datestamps = set()
        if run:
            if run not in self.runs:
                raise KeyError("%s: unrecognised run name" % run)
            runs = [run]
        else:
            runs = self.runs
        for run in runs:
            # For each run, locate the 'flowcell_basecalls.tsv'
            # metadata file
            tsv_file = os.path.join(self.path,
                                    self.run_dirs[run],
                                    "flowcell_basecalls.tsv")
            if os.path.exists(tsv_file):
                # Extract datastamp from each flow cell subdirectory
                flowcell_data = TabFile(tsv_file, first_line_is_header=True)
                for data in flowcell_data:
                    datestamp = get_flow_cell_datestamp(os.path.basename(data['SubDir']))
                    if datestamp:
                        # Only add non-empty values
                        datestamps.add(datestamp)
        try:
            return sorted(list(datestamps))[0]
        except IndexError:
            return None

    @staticmethod
    def _make_project_id(name):
        """
        Internal: generates a project ID

        Arguments:
          name (str): project name to generate ID from
        """
        # Expect name of the form "PromethION_Project_009_PerGynt"
        project_number = RE_PROJECT_DIR_NAME.match(name).group(1)
        return f"PROMETHION#{project_number}"


class ProjectInfo(MetadataDict):
    """
    Class for storing and handling project information

    Defines the following fields:

    - id: project ID
    - name: project name
    - user: user associated with the project
    - PI: PI associated with the project
    - application: application (experiment) type
    - organism: list of one or more organisms
    - data_dir: the source data directory
    - comments: free-text comments field

    Fields are accessed using either e.g.:

    >>> ProjectInfo(filein='project.info').id

    or:

    >>> ProjectInfo(filein='project.info')['id']

    The second form is also used for setting values:

    >>> ProjectInfo(filein='project.info')['user'] = "Per Gynt"

    Data is written to file using the 'save()' method.

    Arguments:
      filein (str): path to exisiting project information
        file to read data from
    """

    def __init__(self, filein=None):
        MetadataDict.__init__(self,
                              attributes={
                                  'name': 'Project name',
                                  'id': 'Project ID',
                                  'platform': 'Platform',
                                  'user': 'User',
                                  'PI': 'PI',
                                  'application': 'Application',
                                  'organism': 'Organism',
                                  'runs': 'Runs',
                                  'data_dir': 'Data directory',
                                  'comments': 'Comments',
                              },
                              order=(
                                  'name',
                                  'id',
                                  'platform',
                                  'user',
                                  'PI',
                                  'application',
                                  'organism',
                                  'runs',
                                  'data_dir',
                                  'comments',
                              ),
                              filen=filein)


class SamplesInfo(MetadataTabFile):
    """
    Class for handling and storing information about samples

    Each sample is stored as an entry with the following
    fields:

    - sample: sample name (must be unique)
    - barcode: barcode ID
    - flowcell: flow cell ID

    Samples are added using the 'add_sample' method, and
    modified using the 'update_sample' method.

    Iterating over a SamplesInfo instance returns each of
    the sample entries:

    >>> for s in SamplesInfo(filein=...):
    ...     print(s)

    Data are written to file using the 'save' method.

    Arguments:
      filein: (optional) name of an existing file to read in
        samples information from
    """

    def __init__(self, filein=None):
        MetadataTabFile.__init__(self,
                                 # Default fields
                                 fields=['Sample',
                                         'Barcode',
                                         'Flowcell'],
                                 # Sort function
                                 sort_func=self._sort_key,
                                 filein=filein)

    def add_sample(self, sample, barcode, flowcell):
        """
        Add an entry for a sample

        Sample must not already appear in the file.
        Wraps the 'add_entry' method of the base
        class.

        Arguments:
          sample (str): name of sample to add
          barcode (str): associated barcode ID
          flowcell (str): associated flowcell ID
        """
        self.add_entry(sample, barcode=barcode, flowcell=flowcell)

    def update_sample(self, sample, **kws):
        """
        Update existing entry for a sample

        Sample must already be present in the file.
        Wraps the 'update_entry' method of the base
        class.

        The 'kws' argument maps data items to the
        new values that should be assigned:

        - barcode: new associated barcode ID
        - flowcell: new associated flowcell ID

        Arguments:
          sample (str): name of sample to update
          kws (mapping): maps data items to new values
        """
        self.update_entry(sample, **kws)

    @staticmethod
    def _sort_key(entry):
        """
        Internal: generates key for sorting entries
        """
        name = entry['Sample']
        return extract_prefix(name), extract_index(name)


class FlowcellBasecallsInfo(TabFile):
    """
    Class for handling and storing info about flow cell basecalls data

    Data about each flow cell basecalls directory is stored as an
    entry with the following fields:

    - run: parent run (maybe 'None')
    - pool_name: name of parent sample pool
    - sub_dir: path to flow cell data directory
    - flow_cell_id: flow cell ID
    - reports: 'True' if there is an HTML report
    - kit: the kit name used to generate the data
    - modifications: modifications used if modified basecalling
      was performed
    - trim_barcodes: whether barcodes were trimmed
    - minknow_version: version of the MinKNOW software used

    Flow cell directories are added using the 'add_base_calls'
    method.

    Iterating over a FlowcellBasecallsInfo instance returns each of
    the basecalls entries:

    >>> for d in FlowcellBasecallsInfo(filein=...):
    ...     print(d)

    Data are written to file using the 'save' method.

    Arguments:
      filein: (optional) name of an existing file to read in
        base calls information from
    """

    def __init__(self, filein=None):
        self._fields = ("Run",
                        "SubDir",
                        "FlowCellID",
                        "Reports",
                        "Kit",
                        "Modifications",
                        "TrimBarcodes",
                        "MinknowVersion",
                        "BasecallingModel",
                        "FileTypes")
        self._kws = tuple([convert_field_name(f)
                           for f in self._fields])
        TabFile.__init__(self,
                         column_names=self._fields,
                         first_line_is_header=True,
                         filen=filein)

    def add_base_calls(self, **kws):
        """
        Add information about a basecalls directory

        Arguments:
          kws (mapping): keywords and values, where keywords
            can be any of the allowed data fields
        """
        for k in kws:
            if k not in self._kws:
                raise KeyError("'%s': unrecognised field" % k)
        data = [kws[k] if k in kws else None for k in self._kws]
        self.append(data)

    def save(self, fileout):
        """
        Save the data to file

        Arguments:
          fileout (str): path to the output file
        """
        self.write(fileout, include_header=True)
