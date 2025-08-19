#!/usr/bin/env python3
#
#     mock.py: utilities for creating mock data for tests
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs
#

"""
Utilities for creating mock data for testing
"""
import json
from pathlib import Path
from textwrap import dedent

from .analysis import ProjectInfo
from .analysis import SamplesInfo
from .analysis import FlowcellBasecallsInfo

# Example JSON data from HTML files

# MinKNOW v24
HTML_JSON_DATA_24 = {
    "run_setup": [{"title": "Flow cell type", "value": "FLO-PRO114M"},
                  {"title": "Flow cell type alias", "value": "FLO-PRO114M"},
                  {"title": "Flow cell ID", "value": "PAW15677"},
                  {"title": "Kit type", "value": "SQK-RBK114-24"}],
    "run_settings": [{"title": "Run limit", "value": "72 hrs"},
                     {"title": "Active channel selection", "value": "On"},
                     {"title": "Pore scan freq.", "value": "1.5 hrs"},
                     {"title": "Reserved pores", "value": "On"},
                     {"title": "Minimum read length", "value": "200 bp"},
                     {"title": "Basecalling", "value": "Super-accurate basecalling, 400 bps"},
                     {"title": "Modified basecalling", "value": "On"},
                     {"title": "Modified base context", "value": "5mC & 5hmC"},
                     {"title": "Trim barcodes", "value": "Off"},
                     {"title": "Mid-read barcode filtering", "value": "Off"},
                     {"title": "Min Q score", "value": "10"}],
    "data_output_settings": [{"title": "FAST5 output", "value": "Off"},
                             {"title": "FASTQ data output", "value": "One file every 10 minutes"},
                             {"title": "POD5 data output", "value": "One file per hour"},
                             {"title": "BAM file output", "value": "On"},
                             {"title": "BAM data output", "value": "One file every 10 minutes"},
                             {"title": "Bulk file output", "value": "Off"},
                             {"title": "Data location", "value": "__DATA_LOCATION__"}],
    "software_versions": [{"title": "MinKNOW", "value": "24.02.19"},
                          {"title": "Bream", "value": "7.9.8"},
                          {"title": "Configuration", "value": "5.9.18"},
                          {"title": "Dorado", "value": "7.3.11"},
                          {"title": "MinKNOW Core", "value": "5.9.12"}],
}

# MinKNOW v25
HTML_JSON_DATA_25 = {
    "run_setup": [{"title": "Flow cell type", "value": "FLO-PRO114M"},
                  {"title": "Flow cell type alias", "value": "FLO-PRO114M"},
                  {"title": "Flow cell ID", "value": "PBC32212"},
                  {"title": "Kit type", "value": "SQK-PCB114-24"}],
    "run_settings": [{"title": "Run limit", "value": "72 hrs"},
                     {"title": "Pore scan freq.", "value": "1.5 hrs"},
                     {"title": "Reserved pores", "value": "On"},
                     {"title": "Basecalling", "value": "High-accuracy model 400bps"},
                     {"title": "Modified basecalling", "value": "Off"},
                     {"title": "Trim barcodes", "value": "Off"},
                     {"title": "Min Q score", "value": "9"}],
    "data_output_settings": [{"title": "FAST5 output", "value": "Off"},
                             {"title": "FASTQ data output", "value": "One file every 10 minutes"},
                             {"title": "POD5 data output", "value": "One file per hour, or 500000000 bases per batch"},
                             {"title": "BAM file output", "value": "On"},
                             {"title": "BAM data output", "value": "One file every 10 minutes"},
                             {"title": "Bulk file output", "value": "Off"},
                             {"title": "Data location", "value": "__DATA_LOCATION__"}],
    "software_versions": [{"title": "MinKNOW", "value": "25.03.7"},
                          {"title": "Bream", "value": "8.4.4"},
                          {"title": "Configuration", "value": "6.4.10"},
                          {"title": "Dorado", "value": "7.8.3"},
                          {"title": "MinKNOW Core", "value": "6.4.8"}],
}


class MockPromethionDataDir:
    """
    Utility class for creating directories which mimic PromethION
    data directories

    Arguments:
        name (str): name of the PromethION project
    """
    def __init__(self, name):
        self.name = str(name)
        self._flow_cells = list()
        self._basecalls_dirs = list()

    def add_flow_cell(self, name, run=None, pool=None):
        """
        Define a flow cell within the mock data

        Arguments:
            name (str): name of the flow cell
            run (str): optional, name of the parent run
            pool (str): optional, name of the parent pool
        """
        self._flow_cells.append(MockFlowcellDir(name, run=run, pool=pool))

    def add_basecalls_dir(self, relpath, flow_cell_name=None):
        """
        Define a basecalls directory within the mock data

        Arguments:
            relpath (str): path to the basecalls dir
              relative to the project directory
            flow_cell_name (str): optional flow cell name
        """
        self._basecalls_dirs.append(MockBasecallsDir(relpath,
                                                     flow_cell_name=flow_cell_name))

    def create(self, top_dir):
        """
        Builds a mock PromethION output directory

        Arguments:
            top_dir (str): path to directory to create the mock
              directory structure under
        """
        # Make top level directory
        top_dir = Path(top_dir).absolute().joinpath(self.name)
        top_dir.mkdir()
        print(f"...made {top_dir}")
        # Add flow cells
        for fc in self._flow_cells:
            fc.create(top_dir)
        # Add basecall dirs
        for bc in self._basecalls_dirs:
            bc.create(top_dir)
        return str(top_dir)


class MockFlowcellDir:
    """
    Utility class for faking flow cell data directories
    """
    def __init__(self, name, run=None, pool=None):
        self.name = str(name)
        self.run = str(run) if run else run
        self.pool = str(pool) if pool else pool

    def create(self, top_dir):
        """
        Create the mock flow cell directory

        Arguments:
            top_dir (Path): directory to create the
              mock flow cell directory under
        """
        fc_path = Path(top_dir)
        if self.run:
            fc_path = fc_path.joinpath(self.run)
        if self.pool:
            fc_path = fc_path.joinpath(self.pool)
        fc_path = fc_path.joinpath(self.name)
        fc_path.mkdir(parents=True)
        print(f"...made {fc_path}")
        # Make subdirs
        self.create_pod5_dirs(fc_path)
        for name in ("bam", "fastq"):
            self.create_data_dirs(name,fc_path)
        self.create_report(fc_path)
        return str(fc_path)

    def create_data_dirs(self, name, top_dir):
        """
        Create data directories (BAM, FASTQ)

        Arguments:
            name (str): name of the data type ("bam", "fastq")
            top_dir (Path): directory to create mock
              data directories under
        """
        for ext in ("pass", "fail"):
            data_dir = top_dir.joinpath(f"{name}_{ext}")
            data_dir.mkdir()
            print(f"...made {data_dir}")
            create_barcode_dirs(data_dir)

    def create_pod5_dirs(self, top_dir):
        """
        Create POD5 directories

        Argument:
          top_dir (Path): directory to create mock
            POD5 directories under
        """
        for name in ("pod5", "pod5_skip"):
            pod5_path = top_dir.joinpath(name)
            pod5_path.mkdir()
            print(f"...made {pod5_path}")

    def create_report(self, top_dir):
        """
        Create HTML report

        Arguments:
          top_dir (Path): directory to create the
            mock report under
        """
        report_file = top_dir.joinpath(f"report_{self.name}.html")
        create_html_report(str(report_file))
        print(f"...made {report_file}")


class MockBasecallsDir:
    """
    Utility class for faking basecalls data directories

    Arguments:
        relpath (str): relative path to the basecalls directory
        flow_cell_name (str): optional, associated flow cell name
    """
    def __init__(self, relpath, flow_cell_name=None):
        self.relpath = str(relpath)
        self.flow_cell_name = flow_cell_name

    def create(self, top_dir):
        """
        Create the mock base calls directory

        Arguments:
            top_dir (Path): directory to create the
              mock flow cell directory under
        """
        bc_path = Path(top_dir).joinpath(self.relpath)
        bc_path.mkdir(parents=True)
        print(f"...made {bc_path}")
        # Make subdirs
        self.create_pass_dir(bc_path)
        self.create_report(bc_path)
        return str(bc_path)

    def create_pass_dir(self, top_dir):
        """
        Create the "pass" subdirectory

        Arguments:
            top_dir (Path): directory to create the
              mock "pass" directory under
        """
        pass_dir = top_dir.joinpath("pass")
        pass_dir.mkdir()
        print(f"...made {pass_dir}")
        create_barcode_dirs(pass_dir)

    def create_report(self, top_dir):
        """
        Create HTML report

        Arguments:
          top_dir (Path): directory to create the
            mock report under
        """
        if self.flow_cell_name:
            report_file = top_dir.joinpath(f"report_{self.flow_cell_name}.html")
            create_html_report(str(report_file))
            print(f"...made {report_file}")


class MockProjectAnalysisDir:
    """
    Utility class for creating directories which mimic PromethION
    project analysis directories

    Arguments:
        name (str): name of the PromethION project analysis
          directory
    """
    def __init__(self, name):
        self.name = str(name)

    def create(self, top_dir, data_dir=None , user=None, principal_investigator=None,
               application=None, organism=None, run_id=None):
        """
        Create mock PromethION analysis project directory

        Arguments:
          top_dir (Path): directory to create the
            mock directory under
          data_dir (str): path to primary data directory
            (optional)
          user (str): optional, name of user
          principal_investigator (str): optional, name of
            principal investigator
          application (str): optional, application
          organism (str): optional, organism
          run_id (str): optional, run ID string
        """
        # Make top level directory
        top_dir = Path(top_dir).absolute().joinpath(self.name)
        top_dir.mkdir()
        print(f"...made {top_dir}")
        # Add in subdirs
        for subdir in ("logs", "ScriptCode"):
            top_dir.joinpath(subdir).mkdir()
        # Add in placeholder files
        for file_name in ("README",):
            top_dir.joinpath(file_name).touch()
        # Add in project info file
        project_info = ProjectInfo()
        metadata_values = {
            'user': user,
            'PI': principal_investigator,
            'application': application,
            'organism': organism,
            'id': run_id
        }
        for k in metadata_values:
            value = metadata_values[k]
            if value:
                project_info[k] = value
        if data_dir:
            project_info['name'] = str(Path(data_dir).name)
            project_info['data_dir'] = str(data_dir)
        project_info['platform'] = "promethion"
        project_info.save(filen=str(top_dir.joinpath("project.info")))
        # Add in samples file
        samples =  SamplesInfo()
        for name, flowcell, barcode in (("PG1", "PAW14589", "NB03"),
                                        ("PG2", "PAW14589", "NB04"),):
            samples.add_sample(name, flowcell, barcode)
        samples.save(fileout=str(top_dir.joinpath("samples.tsv")))
        # Add in flowcells data file
        if data_dir:
            flowcells = FlowcellBasecallsInfo()
            flowcells.save(fileout=str(top_dir.joinpath(f"{str(Path(data_dir).name)}.tsv")))
        return str(top_dir)

def create_barcode_dirs(top_dir, number_of_barcodes=24):
    """
    Creates 'barcode01'...'barcode24' subdirectories

    NB directories are not populated

    Arguments:
      top_dir (Path): directory to create mock
        barcode directories under
      number_of_barcodes (int): number of barcode
        directories to create
    """
    for n in range(0,number_of_barcodes):
        barcode_dir = top_dir.joinpath(f"barcode{n+1:02d}")
        barcode_dir.mkdir()
        # Don't populate with files for now

def create_html_report(file_name, minknow_version="25"):
        """
        Create a mock MinKNOW HTML report file

        Arguments:
          file_name (str): file name and path for new
            mock HTML report file
          minknow_version (str): version of MinKNOW to
            mimick (either "24" or "25"; default is
            "25")
        """
        # Fetch appropriate example JSON data
        if minknow_version == "25":
            json_data = HTML_JSON_DATA_25
        elif minknow_version == "24":
            json_data = HTML_JSON_DATA_24
        else:
            raise Exception(f"Unable to create report for MinKNOW "
                            f"version '{minknow_version}'")
        # Update the "Data location" value
        for s in json_data["data_output_settings"]:
            if s["title"] == "Data location":
                s["value"] == f"\"{str(Path(file_name).parent)}\""
        # Write the mock HTML file
        Path(file_name).write_text(
            dedent(f"""
            const reportDataJson = {json.dumps(json_data)};
            """))
