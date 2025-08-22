#!/usr/bin/env python3
#
#     nanopore.promethion: library for for managing ONT PromethION data
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs
#

"""
Library for exploring & managing Oxford Nanopore Technologies
PromethION data for BCF.

Provides the following classes:

* ProjectDir: handles a PromethION project directory
* FlowCell: handles a PromethION flow cell directory
* BasecallsDir: handles a Promethion basecalls directory
* HtmlReport: handles a MinKNOW HTML report
* JsonReport: handles a MinKNOW JSON report
* BasecallsMetadata: metadata extracted from MinKNOW report
* ProjectAnalysisDir: analysis directory associated with a project

Provides the following functions:

* is_flow_cell_name: checks if string looks a flow cell name
* get_flow_cell_id: extract flow cell ID from a name
* get_flow_cell_datestamp: extract datestamp from a name
* barcode_dirs: get list of "barcode" subdirectories
"""

import os
import re
import json
import shutil
import logging
from pathlib import Path

# Module specific logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RE_FLOW_CELL_NAME = re.compile("^([0-9]+)_([0-9]+)_([0-9][A-Z])_(P[A-Z][A-Z][0-9]+)_([a-z0-9]+)$")


class ProjectDir:
    """
    Class representing a PromethION project directory

    Arguments:
      path (str): path to the project directory
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.name = os.path.basename(self.path)
        # Walk the directory structure looking for "pod5", "*_pass"
        # and "pass" directories
        # Use these to identify flow cell and basecalls directories
        self.flow_cells = []
        self.basecalls_dirs = []
        print("...scanning %s" % self.path)
        flow_cell_dirs = set()
        basecalls_dirs = set()
        for root, dirs, files in os.walk(path):
            for d in dirs:
                if d == "pass":
                    # Base calls directory
                    basecalls_dirs.add(root)
                    break
                elif d == "pod5" or d.endswith("_pass"):
                    # Flow cell directory
                    flow_cell_dirs.add(root)
                    break
        print("...located %d flow cell directories" % len(flow_cell_dirs))
        for d in flow_cell_dirs:
            print("...analysing flow cell '%s'" % d)
            # Match flow cell to parent pool
            fc = FlowCell(d, project_dir=self.path)
            self.flow_cells.append(fc)
        self.flow_cells = sorted(self.flow_cells, key=lambda x: x.name)
        print("...located %d base calls directories" % len(basecalls_dirs))
        for d in basecalls_dirs:
            print("...analysing base call dir '%s'" % d)
            # Look for a matching pool
            for fc in self.flow_cells:
                if fc.pool in d.split(os.sep):
                    bc = BasecallsDir(d, pool=fc.pool, run=fc.run)
                    self.basecalls_dirs.append(bc)
                    break
        self.basecalls_dirs = sorted(self.basecalls_dirs, key=lambda x: x.name)


class FlowCell:
    """
    Class representing a PromethION flow cell directory

    Arguments:
      path (str): path to the flow cell directory
      project_dir (str): optional, path to the parent
        PromethION project directory
    """

    def __init__(self, path, project_dir=None):
        # Absolute path
        self.path = os.path.abspath(path)
        self.name = os.path.basename(self.path)
        # Assign parent pool and run
        if not project_dir:
            self.project_dir = None
        else:
            self.project_dir = os.path.abspath(project_dir)
        pool_dir = os.path.dirname(self.path)
        self.pool = os.path.basename(pool_dir)
        run_dir = os.path.dirname(pool_dir)
        # Check parent pool or run aren't same as project dir
        if self.project_dir is not None and self.project_dir in (pool_dir, run_dir):
            logger.warning("%s: missing run or pool directories?" % self.path)
            self.run = None
        else:
            self.run = os.path.basename(run_dir)
        # Associated metadata
        self.metadata = BasecallsMetadata()
        # Name and ID for flow cell
        if not is_flow_cell_name(self.name):
            raise Exception("'%s': not a flow cell name?" % self.name)
        self.id = get_flow_cell_id(self.name)
        self.datestamp = get_flow_cell_datestamp(self.name)
        # POD5 directory
        for d in ("pod5", "pod5_pass"):
            self.pod5 = os.path.join(self.path, d)
            if os.path.exists(self.pod5):
                break
            else:
                self.pod5 = None
        # BAMs directory
        self.bam_pass = os.path.join(self.path, "bam_pass")
        if not os.path.exists(self.bam_pass):
            self.bam_pass = None
        # FASTQs directory
        self.fastq_pass = os.path.join(self.path, "fastq_pass")
        if not os.path.exists(self.fastq_pass):
            self.fastq_pass = None
        # Reports
        self.reports = []
        self.sample_sheet = None
        for f in os.listdir(self.path):
            if f.startswith("report_"):
                self.reports.append(f)
                if f.endswith(".html") or f.endswith(".json"):
                    try:
                        self.metadata.load_from_report(
                            os.path.join(self.path, f))
                    except Exception as ex:
                        print(f"{f}: failed to load metadata from file "
                              f"(ignored): {ex}")
            elif f.startswith("sample_sheet_"):
                self.sample_sheet = f

    @property
    def file_types(self):
        file_types = []
        if self.pod5:
            file_types.append("pod5")
        if self.bam_pass:
            file_types.append("bam")
        if self.fastq_pass:
            file_types.append("fastq")
        return file_types

    @property
    def report_types(self):
        report_types = []
        if self.html_report:
            report_types.append("html")
        if self.json_report:
            report_types.append("json")
        return report_types

    @property
    def html_report(self):
        for f in self.reports:
            if f.endswith(".html"):
                return os.path.join(self.path, f)
        return None

    @property
    def json_report(self):
        for f in self.reports:
            if f.endswith(".json"):
                return os.path.join(self.path, f)
        return None

    def __repr__(self):
        return "%s/%s" % (self.pool, self.name)


class BasecallsDir:
    """
    Class representing a PromethION basecalls directory

    Arguments:
      path (str): path to the basecalls directory
      pool (str): optional, pool name to associate the
        basecalls directory with
      run (str): optional, run name to associate the
        basecalls directory
    """

    def __init__(self, path, pool=None, run=None):
        self.path = os.path.abspath(path)
        self.name = os.path.basename(path)
        self.parent = os.path.basename(os.path.dirname(self.path))
        # Assign parent pool and run
        self.pool = pool
        self.run = run
        # Associated metadata
        self.metadata = BasecallsMetadata()
        # BAMs and FASTQs directory
        self.pass_dir = os.path.join(self.path, "pass")
        if not os.path.exists(self.pass_dir):
            self.pass_dir = None
        # Reports
        self.reports = []
        self.sample_sheet = None
        for f in os.listdir(self.path):
            if f.startswith("report_"):
                self.reports.append(f)
                if f.endswith(".html") or f.endswith(".json"):
                    try:
                        self.metadata.load_from_report(
                            os.path.join(self.path, f))
                    except Exception as ex:
                        print(f"{f}: failed to load metadata from file "
                              f"(ignored): {ex}")
            elif f.startswith("sample_sheet_"):
                self.sample_sheet = f

    @property
    def file_types(self):
        if self.pass_dir:
            return ["bam", "fastq"]
        else:
            return []

    @property
    def report_types(self):
        report_types = []
        if self.html_report:
            report_types.append("html")
        if self.json_report:
            report_types.append("json")
        return report_types

    @property
    def html_report(self):
        for f in self.reports:
            if f.endswith(".html"):
                return os.path.join(self.path, f)
        return None

    @property
    def json_report(self):
        for f in self.reports:
            if f.endswith(".json"):
                return os.path.join(self.path, f)
        return None

    def __repr__(self):
        return "%s/%s" % (self.parent, self.name)


class HtmlReport:
    """
    Class for handling HTML report from MinKNOW

    Arguments:
      path (str): path to the HTML file
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def extract_json(self):
        """
        Extract and return JSON data embedded in report

        Attempts to locate the JSON data embedded in an HTML
        report file and return as a JSON object.
        """
        json_data = None
        with open(self.path, "rt") as fp:
            for line in fp:
                if line.strip().startswith("const reportDataJson = {"):
                    json_data = line.strip()[len("const reportDataJson = "):-1]
                    break
                elif line.strip().startswith("const reportData={"):
                    json_data = line.strip()[len("const reportData="):]
                    break
        if json_data is None:
            raise Exception("%s: unable to extract JSON data" % self.path)
        return json.loads(json_data)

    def __repr__(self):
        return self.path


class JsonReport:
    """
    Class for handling JSON report from MinKNOW

    Arguments:
      path (str): path to the JSON file
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def extract_json(self):
        """
        Read and return JSON data from the file

        Returns data as a JSON object.
        """
        with open(self.path, "rt") as fp:
            try:
                return json.load(fp)
            except Exception as ex:
                raise Exception(f"{self.path}: unable to extract "
                                f"JSON data: {ex}")

    def __repr__(self):
        return self.path


class BasecallsMetadata:
    """
    Class representing data about a set of basecalls

    The data is populated by extracting data from an HTML
    and JSON report files (by invoking the 'load_from_report'
    method on each file).
    """

    def __init__(self):
        self.html_json_data = None
        self.json_data = None
        self.flow_cell_id = None
        self.flow_cell_type = None
        self.kit = None
        self.basecalling = None
        self.modified_basecalling = None
        self.modifications = None
        self.trim_barcodes = None
        self.software_versions = None
        self.basecalling_model = None
        self.basecalling_config = None

    def load_from_report(self, report_file):
        """
        Load metadata from a report file
        """
        if report_file.endswith(".html"):
            return self.load_from_report_html(report_file)
        elif report_file.endswith(".json"):
            return self.load_from_report_json(report_file)
        else:
            raise Exception(f"{report_file}: unsupported report file "
                            f"type")

    def load_from_report_html(self, html_file):
        """
        Load metadata from an HTML report file

        Arguments:
          html_file (str): path to the HTML report
        """
        self.html_json_data = HtmlReport(html_file).extract_json()
        data = {}
        for k in ('run_settings',
                  'run_setup',
                  'software_versions'):
            data[k] = self._extract_section(self.html_json_data, k)
        # Set values
        setup = data['run_setup']
        ##print(setup)
        self.flow_cell_id = setup['flow_cell_id']
        self.flow_cell_type = setup['flow_cell_type']
        self.kit = setup['kit_type']
        settings = data['run_settings']
        ##print(settings)
        self.basecalling = settings['basecalling']
        self.modified_basecalling = settings['modified_basecalling']
        if self.modified_basecalling == "On":
            for name in ("modifications", "modified_base_context"):
                self.modifications = self._fetch_item(settings, name)
                if self.modifications:
                    break                         
        else:
            self.modifications = None
        try:
            self.trim_barcodes = settings['trim_barcodes']
        except KeyError:
            self.trim_barcodes = None
        # Store the software versions
        self.software_versions = { k:data["software_versions"][k]
                                   for k in data["software_versions"] }
        # Finish
        return self

    def html_json(self):
        """
        Return JSON data extracted from the HTML report
        """
        if self.html_json_data is None:
            return {}
        return json.dumps(self.html_json_data, sort_keys=True,
                          indent=4)

    def load_from_report_json(self, json_file):
        """
        Load metadata from a JSON report file

        Arguments:
          json_file (str): path to the JSON report
        """
        self.json_data = JsonReport(json_file).extract_json()
        try:
            for acquisition in self.json_data["acquisitions"]:
                # Basecalling model
                if not self.basecalling_model:
                    try:
                        self.basecalling_model = acquisition["acquisition_run_info"]["config_summary"]["basecalling_model_version"]
                    except KeyError:
                        pass
                # Basecalling config
                if not self.basecalling_config:
                    try:
                        self.basecalling_config = acquisition["acquisition_run_info"]["config_summary"]["basecalling_config_filename"]
                    except KeyError:
                        pass
        except KeyError:
            pass
        return self

    def json(self):
        """
        Return data extracted from the JSON report
        """
        if self.json_data is None:
            return {}
        return json.dumps(self.json_data, sort_keys=True,
                          indent=4)

    def _extract_section(self, json_data, name):
        """
        Internal: extracts data from JSON for a section

        Given a section name in the JSON data, returns a
        dictionary where the 'title' of each item is a
        key with the corresponding 'value'.

        Key names are generated from 'title' strings by
        converting to lower case and replacing ' ', '-'
        and '.' characters to underscores.

        Arguments:
          json_data (dict): JSON data
          name (str): name of the section to extract
        """
        return {
            s['title'].lower()
            .replace(' ', '_')
            .replace('-', '_')
            .replace('.', ''): s['value']
            for s in json_data[name]
        }

    def _fetch_item(self, data, name):
        """
        """
        try:
            return data[name]
        except KeyError:
            logger.warning(f"'{name}': metadata item not found")
            return None


def is_flow_cell_name(name):
    """
    Check if a name looks like a flow cell

    Arguments:
      name (str): putative flow cell name

    Returns:
      Boolean: True if name is a flow cell name,
        False otherwise.
    """
    return bool(RE_FLOW_CELL_NAME.match(name))


def get_flow_cell_id(name):
    """
    Extract the flow cell ID from a name

    Arguments:
      name (str): flow cell name

    Return:
      String: flow cell ID, or "" if ID cannot be extracted.
    """
    try:
        return RE_FLOW_CELL_NAME.match(name).group(4)
    except AttributeError:
        return ""


def get_flow_cell_datestamp(name):
    """
    Extract the flow cell date stamp from a name

    Arguments:
      name (str): flow cell name

    Return:
      String: flow cell date stamp.
    """
    try:
        return RE_FLOW_CELL_NAME.match(name).group(1)
    except AttributeError:
        return ""


def barcode_dirs(d):
    """
    Return list of 'barcode' subdirectories

    Arguments:
      d (str): path to the directory to examine

    Returns:
      List: list of subdirectory names under the supplied
        directory with names 'barcode...'.
    """
    return sorted([b for b in filter(lambda x: x.startswith("barcode")
                   and os.path.isdir(os.path.join(d, x)),
                                     os.listdir(d))])
