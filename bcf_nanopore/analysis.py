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
from bcftbx.TabFile import TabFile
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from .nanopore.promethion import ProjectDir
from .utils import MetadataTabFile
from .utils import convert_field_name
from .utils import fmt_value
from .utils import fmt_yes_no

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RE_PROJECT_DIR_NAME = re.compile("^PromethION_Project_([0-9]+)_(.+)$")


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
          project.info
          samples.tsv
    
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
        # Sample sheet (per sample metadata)
        self.samples_file = os.path.join(self.path, "samples.tsv")
        if os.path.exists(self.samples_file):
            self.samples_info = SamplesInfo(self.samples_file)
        else:
            logger.warning("%s: no 'samples.tsv' file found" %
                           self.path)
            self.samples_info = SamplesInfo()

    def create(self, project_dir, user, PI, application, organism,
               sample_sheet=None, samples=[]):
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
          samples (list): optional, list of tuples
            of the form (SAMPLE, BARCODE, FLOWCELL)
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
        # Samples
        # Copy in data from external sample sheet
        # Sample sheet should be a CSV format file with an
        # initial header line (which is ignored) followed by
        # lines with either 2 or 3 fields
        if samples:
            for sample in samples:
                sample_name, barcode, flowcell = sample
                self.samples_info.add_sample(sample_name,
                                             barcode,
                                             flowcell)
            self.samples_info.save(fileout=self.samples_file)
        # Collate flow cell and base calls information into TSV
        project = ProjectDir(project_dir)
        flow_cells_file = os.path.join(self.path,
                                       "%s.tsv" % self.info.name)
        fc_file = FlowcellBasecallsInfo()
        for fc in project.flow_cells:
            fc_file.add_base_calls(
                run=fc.run,
                pool_name=fc.pool,
                sub_dir=fc,
                flow_cell_id=fc.id,
                reports=fmt_yes_no(fc.html_report),
                kit=fmt_value(fc.metadata.kit),
                modifications=("none"
                               if fc.metadata.modified_basecalling == "Off"
                               else fmt_value(fc.metadata.modifications)),
                trim_barcodes=fmt_value(fc.metadata.trim_barcodes))
        for bc in project.basecalls_dirs:
            fc_file.add_base_calls(
                run=bc.run,
                pool_name=bc.name,
                sub_dir=bc,
                flow_cell_id=fmt_value(bc.metadata.flow_cell_id),
                reports=fmt_yes_no(bc.html_report),
                kit=fmt_value(bc.metadata.kit),
                modifications=("none"
                               if bc.metadata.modified_basecalling == "Off"
                               else fmt_value(bc.metadata.modifications)),
                trim_barcodes=fmt_value(bc.metadata.trim_barcodes))
        fc_file.save(flow_cells_file)
        # Get the earliest date stamp from flow cell names
        try:
            datestamps = set()
            for fc in project.flow_cells:
                datestamps.add(fc.datestamp)
            datestamps = sorted(list(datestamps))
            self.info['datestamp'] = datestamps[0]
            self.info.save(filen=self.project_info_file)
        except Exception as ex:
            print(f"Failed to set datestamp: {ex}")
        # Make subdirectories
        for subdir in ("logs", "ScriptCode"):
            Path(self.path).joinpath(subdir).mkdir()
        # Copy in reports
        reports_dir = os.path.join(self.path, "reports")
        os.makedirs(reports_dir)
        for fc in project.flow_cells:
            report = fc.html_report
            if not report:
                continue
            target = os.path.join(reports_dir,
                                  "%s_%s_%s" % (fc.pool,
                                                fc.id,
                                                os.path.basename(report)))
            shutil.copy(report, target)
        for bc in project.basecalls_dirs:
            report = bc.html_report
            if not report:
                continue
            target = os.path.join(reports_dir,
                                  "%s_%s_%s_%s" % (bc.parent,
                                                   bc.pool,
                                                   bc.metadata.flow_cell_id,
                                                   os.path.basename(report)))
            shutil.copy(report, target)
        # Create a README file
        read_me_file = os.path.join(self.path, "README")
        with open(read_me_file, 'wt') as read_me:
            read_me.write(
                f"""This is the analysis directory for {os.path.basename(self.info.data_dir)}

The following files have been automatically generated:

- '{os.path.basename(self.project_info_file)}': top-level information about the project
- '{os.path.basename(flow_cells_file)}': TSV file with information about
  the flow cell and base calling directories (extracted from the primary
  data directory)
""")
            if os.path.exists(self.samples_file):
                read_me.write(
                    f"- '{os.path.basename(self.samples_file)}': TSV file matching sample names "
                    "to flow cell and barcode IDs\n")
            read_me.write(f"""
The 'reports' subdirectory contains copies of the HTML reports that were
found in the primary data directory (renamed to identify the associated
locations).""")
        # Report information
        print(self.report(mode="summary",
                          fields="name,id,datestamp,platform,analysis_dir,,"
                                 "user,pi,application,organism,primary_data,"
                                 "comments"))

    def exists(self):
        """
        Check whether project analysis directory exists
        """
        return os.path.exists(self.path)

    def report(self, mode, fields):
        """
        Report information on the project analysis directory

        Arguments:
            mode (str): reporting mode
            fields (str): comma separated list of field names
        """
        delimiter = {
            'summary': '\n',
            'tsv': '\t'
        }
        if mode not in ("summary", "tsv"):
            raise Exception("%s: unknown reporting mode" % mode)
        # Collect data to report
        fields = fields.split(',')
        names = []
        values = []
        for f in fields:
            field = f.strip().lower()
            fmt_func = fmt_value
            if field == "" or field == "null":
                name = None
                value = ''
            elif field == "name":
                name = "Project name"
                value = self.info.name
            elif field == "id":
                name = "Project ID"
                value = self.info.id
            elif field == "datestamp":
                name = "Datestamp"
                value = self.info.datestamp
            elif field == "platform":
                name = "Platform"
                value = self.info.platform
            elif field == "user":
                name = "User"
                value = self.info.user
            elif field == "pi":
                name = "PI"
                value = self.info.PI
            elif field == "application":
                name = "Application"
                value = self.info.application
            elif field == "organism":
                name = "Organism"
                value = self.info.organism
            elif field == "nsamples" or field == "#samples":
                name = "#samples"
                value = len(self.samples_info)
                def fmt_func(n): return '?' if n == 0 else str(n)
            elif field == 'sample_names' or field == 'samples':
                name = "samples"
                for s in self.samples_info:
                    print(s)
                value = ",".join([s["Sample"] for s in self.samples_info])
                def fmt_func(s): return '?' if s == "" else s
            elif field == "primary_data":
                name = "Primary data dir"
                value = self.info.data_dir
            elif field == "analysis_dir":
                name = "Analysis dir"
                value = self.path
            elif field == "comments":
                name = "Comments"
                value = self.info.comments
                def fmt_func(s): return '' if s is None else str(s)
            else:
                raise KeyError("%s: unrecognised field" % f)
            names.append(name)
            values.append(fmt_func(value))
        # Generate output
        output = []
        if mode == "summary":
            # Add title in summary mode
            output.extend([self.info.name, '=' * len(self.info.name)])
        for name, value in zip(names, values):
            if mode == "summary":
                if name:
                    output.append("%-16s: %s" % (name, value))
                else:
                    output.append('')
            elif mode == "tsv":
                output.append(str(value))
        return delimiter[mode].join(output)

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
                                  'datestamp': 'Datestamp',
                                  'platform': 'Platform',
                                  'user': 'User',
                                  'PI': 'PI',
                                  'application': 'Application',
                                  'organism': 'Organism',
                                  'data_dir': 'Data directory',
                                  'comments': 'Comments',
                              },
                              order=(
                                  'name',
                                  'id',
                                  'datestamp',
                                  'platform',
                                  'user',
                                  'PI',
                                  'application',
                                  'organism',
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
                        "PoolName",
                        "SubDir",
                        "FlowCellID",
                        "Reports",
                        "Kit",
                        "Modifications",
                        "TrimBarcodes")
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
