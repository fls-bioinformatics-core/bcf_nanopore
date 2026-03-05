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

JSON_DATA_24 = {
    "acquisitions": [
        { "acquisition_run_info":
          { "config_summary":
            { "basecalling_config_filename": "dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg",
              "basecalling_enabled": "true",
              "basecalling_model_version": "dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
              "channel_count": 3000 }
           }
         }
    ]
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

JSON_DATA_25 = {
    "acquisitions": [
        { "acquisition_run_info":
          { "config_summary":
            { "barcoding_enabled": "true",
              "barcoding_kits": [ "SQK-NBD114-24" ],
              "basecalling_enabled": "true",
              "basecalling_model_names": {
                  "modified_models": [ "dna_r10.4.1_e8.2_400bps_hac@v4.3.0_5mCG_5hmCG@v1" ],
                  "simplex_model": "dna_r10.4.1_e8.2_400bps_hac@v4.3.0" },
              "basecalling_model_version": "dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
              "channel_count": 3000,
             }
           }
         }
    ]
}


class MockPromethionDataDir:
    """
    Utility class for creating directories which mimic PromethION
    data directories

    Example of creating an initial mock project directory with two
    flow cells

    >>> project = MockPromethionDataDir("PromethION_Project_011_JohnDoe")
    >>> project.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
    ...                       relpath=Path("PG1-2_20240513").joinpath("PG1-2"))
    >>> project.add_flow_cell("20240529_0830_1A_PAW17328_523ce32d",
    ...                       relpath=Path("PG3-4_20240529").joinpath("PG3-4"))
    >>> project_dir = project.create(os.getcwd())

    To add another flow cell to this project directory:

    >>> project.add_flow_cell("20240502_0831_1A_PAW16072_253da29d"
    ...                       relpath=Path("PG5-6_20240602").joinpath("PG5-6"))
    >>> project_dir = project.update(os.getcwd())

    Arguments:
        name (str): name of the PromethION project
    """
    def __init__(self, name):
        self.name = str(name)
        self._flow_cells = list()
        self._basecalls_dirs = list()

    def add_flow_cell(self, name, relpath=None):
        """
        Define a flow cell within the mock data

        Arguments:
            name (str): name of the flow cell
            relpath (str): optional, path to the flow cell
              parent directory relative to the project
              directory
        """
        self._flow_cells.append(MockFlowcellDir(name, relpath=relpath))

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
        mock_dir = Path(top_dir).absolute().joinpath(self.name)
        mock_dir.mkdir()
        print(f"...made {mock_dir}")
        return self.update(top_dir)

    def update(self, top_dir):
        """
        Update a mock PromethION output directory

        Args:
            top_dir (str): path to directory to update the mock
            directory structure under
        """
        # Check top level directory
        mock_dir = Path(top_dir).absolute().joinpath(self.name)
        if not mock_dir.exists():
            raise OSError(f"Directory {mock_dir} does not exist")
        print(f"...updating {mock_dir}")
        # Add flow cells
        for fc in self._flow_cells:
            if not fc.path(mock_dir).exists():
                fc.create(mock_dir)
        # Add basecall dirs
        for bc in self._basecalls_dirs:
            if not bc.path(mock_dir).exists():
                bc.create(mock_dir)
        return str(mock_dir)

class MockFlowcellDir:
    """
    Utility class for faking flow cell data directories

    Arguments:
        name (str): name of the flow cell
        relpath (str): relative path to the flow cell directory
            from the top level directory
    """
    def __init__(self, name, relpath=None):
        self.name = str(name)
        self.relpath = relpath

    def path(self, top_dir):
        """
        Return path to mock flow cell directory

        Args:
            top_dir (Path): directory to create
              mock flow cell directory under

        Returns:
            Path: path to the mock flow cell directory
        """
        if self.relpath:
            return Path(top_dir).joinpath(self.relpath,
                                          self.name)
        return Path(top_dir).joinpath(self.name)

    def create(self, top_dir):
        """
        Create the mock flow cell directory

        Arguments:
            top_dir (Path): directory to create the
              mock flow cell directory under
        """
        fc_path = self.path(top_dir)
        fc_path.mkdir(parents=True)
        print(f"...made {fc_path}")
        # Make subdirs
        self.create_pod5_dirs(fc_path)
        for name in ("bam", "fastq"):
            self.create_data_dirs(name,fc_path)
        self.create_reports(fc_path)
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

    def create_reports(self, top_dir, reports=("html", "json"),
                       minknow_version="25"):
        """
        Create report files

        Arguments:
          top_dir (Path): directory to create the
            mock reports under
          reports (list): list of the report types to create
            (default: "html" and "json")
          minknow_version (str): version of MinKNOW to
            mimick when writing reports (either "24" or
            "25", default: "25")
        """
        for report in reports:
            self.create_report(top_dir, report, minknow_version)

    def create_report(self, top_dir, report_type, minknow_version):
        """
        Create report file

        Arguments:
          top_dir (Path): directory to create the
            mock report under
          report_type (str): either "html" or "json"
          minknow_version (str): either "24" or "25"
        """
        report_file = top_dir.joinpath(f"report_{self.name}.{report_type}")
        if report_type == "html":
            create_html_report(str(report_file), minknow_version=minknow_version)
        elif report_type == "json":
            create_json_report(str(report_file), minknow_version=minknow_version)
        else:
            raise Exception(f"{report_type}: unsupported report type")
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

    def path(self, top_dir):
        """
        Return path to mock basecalls directory

        Args:
            top_dir (Path): directory to create
              mock basecalls directory under

        Returns:
            Path: path to the mock basecalls directory
        """
        return Path(top_dir).joinpath(self.relpath)

    def create(self, top_dir):
        """
        Create the mock base calls directory

        Arguments:
            top_dir (Path): directory to create the
              mock flow cell directory under
        """
        bc_path = self.path(top_dir)
        bc_path.mkdir(parents=True)
        print(f"...made {bc_path}")
        # Make subdirs
        self.create_pass_dir(bc_path)
        self.create_reports(bc_path)
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

    def create_reports(self, top_dir, reports=("html", "json"),
                       minknow_version="25"):
        """
        Create report files

        Arguments:
          top_dir (Path): directory to create the
            mock reports under
          reports (list): list of the report types to create
            (default: "html" and "json")
          minknow_version (str): version of MinKNOW to
            mimick when writing reports (either "24" or
            "25", default: "25")
        """
        for report in reports:
            self.create_report(top_dir, report, minknow_version)

    def create_report(self, top_dir, report_type, minknow_version):
        """
        Create report file

        Arguments:
          top_dir (Path): directory to create the
            mock report under
          report_type (str): either "html" or "json"
          minknow_version (str): either "24" or "25"
        """

        if self.flow_cell_name:
            report_file = top_dir.joinpath(
                f"report_{self.flow_cell_name}.{report_type}")
            if report_type == "html":
                create_html_report(str(report_file),
                                   minknow_version=minknow_version)
            elif report_type == "json":
                create_json_report(str(report_file),
                                   minknow_version=minknow_version)
            else:
                raise Exception(f"{report_type}: unsupported report type")
            print(f"...made {report_file}")


class MockProjectAnalysisDir:
    """
    Utility class for creating directories which mimic PromethION
    project analysis directories

    Create a MockPromethionAnalysisDir instance, add runs to it
    using the add_run() method, then create the mock directory
    structure using the create() method.

    Example usage:

    >>> project = MockProjectAnalysisDir("MyProject")
    >>> project.add_run("Run1", samples={"Sample1": ("barcode01","FC12345")})
    >>> project_dir = project.create("/path/to/top/dir", data_dir="/path/to/data/dir",
    ...                             user="jdoe", principal_investigator="J. Doe",
    ...                             application="Methylation study", organism="E. coli",
    ...                             project_id="PROMETHION#001")

    Arguments:
        name (str): name of the PromethION project analysis
          directory
    """
    def __init__(self, name):
        self.name = str(name)
        self.runs = {}
        self.run_metadata = {}

    def add_run(self, run, samples=None, metadata=None):
        """
        Add a mock run to the analysis project

        If supplied then 'samples' is a dictionary mapping
        sample names to (barcode,flowcell) tuples, for
        example:

        ::

            samples = {
                "Sample1": ("barcode01","FC12345"),
                "Sample2": ("barcode02","FC12345")
            }

        Metadata values can be set by supplying a
        dictionary via the 'metadata' argument which
        maps metadata items to their values, for
        example:

        ::

            metadata = {
                "Order numbers": "#00123",
            }

        These are then written to the 'run.info' file
        in the mock analysis directory.

        Arguments:
          run (str): name of the run to add
          samples (dict): optional, dictionary mapping
            sample names to (barcode,flowcell) tuples
          metadata (dict): optional, dictionary mapping
            items to value for the 'run.info' file
        """
        self.runs[run] = samples
        self.run_metadata[run] = metadata

    def create(self, top_dir, data_dir=None , user=None, principal_investigator=None,
               application=None, organism=None, project_id=None, extra_project_metadata=None):
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
          project_id (str): optional, project ID string
          extra_project_metadata (dict): optional,
            dictionary of extra metadata items and values
            to add to the project metadata
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
            'id': project_id,
            'runs': ",".join(list(self.runs.keys()))
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
        if extra_project_metadata:
            # Append extra items to the metadata file
            with open(str(top_dir.joinpath("project.info")), "at") as fp:
                for key in extra_project_metadata:
                    value = extra_project_metadata[key]
                    if value is None:
                        value = "."
                    fp.write(f"{key}\t{value}\n")
        # Add in run subdirectories
        for ix, run in enumerate(self.runs.keys()):
            run_dir = top_dir.joinpath(f"{ix+1:03d}_{run}")
            run_dir.mkdir()
            print(f"...made run '{run}' ({run_dir})")
            # Add artefacts to run directory
            # Placeholder README file
            with open(run_dir.joinpath("README"), "wt") as fp:
                fp.write("Placeholder README file\n")
            # flowcell_basecalls.tsv file
            with open(run_dir.joinpath("flowcell_basecalls.tsv"), "wt") as fp:
                fp.write("#%s\n" % "\t".join(["Run", "SubDir", "FlowCellID",
                                              "Reports", "Kit", "Modifications",
                                              "TrimBarcodes", "MinknowVersion",
                                              "BasecallingModel", "FileTypes"]))
                if self.runs[run]:
                    flowcells = set()
                    for sample in self.runs[run]:
                        barcode, flowcell = self.runs[run][sample]
                        if flowcell in flowcells:
                            continue
                        flowcells.add(flowcell)
                        subdir = f"20240513_0716_1F_{flowcell}_30105f28"
                        fp.write("%s\n" % "\t".join([run, subdir, flowcell,
                                                     "html", "SQK-PCB114-24", "none",
                                                     "Off", "25.03.7",
                                                     "dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
                                                     "pod5,bam,fastq"]))
            # Samples file
            with open(run_dir.joinpath("samples.tsv"), "wt") as fp:
                fp.write("#Sample\tBarcode\tFlowcell\n")
                if self.runs[run]:
                    for sample in self.runs[run]:
                        barcode, flowcell = self.runs[run][sample]
                        fp.write(f"{sample}\t{barcode}\t{flowcell}\n")
            with open(run_dir.joinpath("samples.tsv"), "rt") as fp:
                print(fp.read())
            # Metadata file
            with open(run_dir.joinpath("run.info"), "wt") as fp:
                fp.write(f"Run name\t{run}\n")
                if self.run_metadata[run]:
                    for item, value in self.run_metadata[run].items():
                        if value is None:
                            value = "."
                        fp.write(f"{item}\t{value}\n")
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

def create_json_report(file_name, minknow_version="25"):
        """
        Create a mock MinKNOW JSON report file

        Arguments:
          file_name (str): file name and path for new
            mock JSON report file
          minknow_version (str): version of MinKNOW to
            mimick (either "24" or "25"; default is
            "25")
        """
        # Fetch appropriate example JSON data
        if minknow_version == "25":
            json_data = JSON_DATA_25
        elif minknow_version == "24":
            json_data = JSON_DATA_24
        else:
            raise Exception(f"Unable to create report for MinKNOW "
                            f"version '{minknow_version}'")
        # Write the mock JSON file
        Path(file_name).write_text(json.dumps(json_data) + "\n")
