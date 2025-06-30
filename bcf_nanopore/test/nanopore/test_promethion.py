#!/usr/bin/env python3

# Tests for the 'nanopore.promethion' module

import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.mock import MockPromethionDataDir
from bcf_nanopore.mock import MockFlowcellDir
from bcf_nanopore.mock import MockBasecallsDir
from bcf_nanopore.mock import create_html_report
from bcf_nanopore.nanopore.promethion import is_flow_cell_name
from bcf_nanopore.nanopore.promethion import get_flow_cell_id
from bcf_nanopore.nanopore.promethion import get_flow_cell_datestamp
from bcf_nanopore.nanopore.promethion import ProjectDir
from bcf_nanopore.nanopore.promethion import FlowCell
from bcf_nanopore.nanopore.promethion import BasecallsDir
from bcf_nanopore.nanopore.promethion import HtmlReport
from bcf_nanopore.nanopore.promethion import BasecallsMetadata


class TestProjectDir(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_project_dir(self):
        """
        ProjectDir: load from directory
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_flow_cell("20240513_0829_1B_PAW15451_20b55ef2", run="PG1-4_20240513", pool="PG3-4")
        project_dir = data_dir.create(self.wd)
        project = ProjectDir(project_dir)
        self.assertEqual(project.name, "PromethION_Project_001_PerGynt")
        self.assertEqual(len(project.flow_cells), 2)

class TestFlowCell(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_flow_cell(self):
        """
        FlowCell: load directly from flow cell directory
        """
        data_dir = MockFlowcellDir("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        flow_cell_dir = data_dir.create(Path(self.wd))
        flow_cell = FlowCell(flow_cell_dir)
        self.assertEqual(flow_cell.name, "20240513_0829_1A_PAW15419_465bb23f")
        self.assertEqual(flow_cell.path, flow_cell_dir)
        self.assertEqual(flow_cell.id, "PAW15419")
        self.assertEqual(flow_cell.datestamp, "20240513")
        self.assertEqual(flow_cell.project_dir, None)
        self.assertEqual(flow_cell.run, "PG1-4_20240513")
        self.assertEqual(flow_cell.pool, "PG1-2")
        self.assertEqual(flow_cell.html_report,
                         str(Path(flow_cell_dir).joinpath("report_20240513_0829_1A_PAW15419_465bb23f.html")))
        self.assertEqual(flow_cell.pod5,str(Path(flow_cell_dir).joinpath("pod5")))
        self.assertEqual(flow_cell.bam_pass, str(Path(flow_cell_dir).joinpath("bam_pass")))
        self.assertEqual(str(flow_cell), "PG1-2/20240513_0829_1A_PAW15419_465bb23f")

    def test_flow_cell_from_project(self):
        """
        FlowCell: load as part of project
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        project_dir = data_dir.create(self.wd)
        flow_cell_dir = Path(project_dir).joinpath("PG1-4_20240513",
                                                    "PG1-2",
                                                    "20240513_0829_1A_PAW15419_465bb23f")
        project = ProjectDir(project_dir)
        flow_cell = project.flow_cells[0]
        self.assertEqual(flow_cell.name, "20240513_0829_1A_PAW15419_465bb23f")
        self.assertEqual(flow_cell.path, str(flow_cell_dir))
        self.assertEqual(flow_cell.id, "PAW15419")
        self.assertEqual(flow_cell.datestamp, "20240513")
        self.assertEqual(flow_cell.project_dir, project_dir)
        self.assertEqual(flow_cell.run, "PG1-4_20240513")
        self.assertEqual(flow_cell.pool, "PG1-2")
        self.assertEqual(flow_cell.html_report,
                         str(flow_cell_dir.joinpath("report_20240513_0829_1A_PAW15419_465bb23f.html")))
        self.assertEqual(flow_cell.pod5,str(flow_cell_dir.joinpath("pod5")))
        self.assertEqual(flow_cell.bam_pass, str(flow_cell_dir.joinpath("bam_pass")))
        self.assertEqual(str(flow_cell), "PG1-2/20240513_0829_1A_PAW15419_465bb23f")

class TestBasecallsDir(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_basecalls_dir(self):
        """
        BasecallsDir: load directly from basecalls directory
        """
        data_dir = MockBasecallsDir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        basecalls_dir = data_dir.create(Path(self.wd))
        basecalls = BasecallsDir(basecalls_dir)
        self.assertEqual(basecalls.path, basecalls_dir)
        self.assertEqual(basecalls.name, "PG1-2")
        self.assertEqual(basecalls.parent, "Rebasecalling")
        self.assertEqual(basecalls.pool, None)
        self.assertEqual(basecalls.run, None)
        self.assertEqual(basecalls.bam_pass, str(Path(basecalls_dir).joinpath("pass")))
        self.assertEqual(basecalls.html_report,
                         str(Path(basecalls_dir).joinpath("report_20240513_0829_1A_PAW15419_465bb23f.html")))

    def test_basecalls_dir_from_project(self):
        """
        BasecallsDir: load directly from basecalls directory
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        project = ProjectDir(project_dir)
        basecalls_dir = str(Path(project_dir).joinpath("PG1-4_20240513","Rebasecalling","PG1-2"))
        basecalls = project.basecalls_dirs[0]
        self.assertEqual(basecalls.path, basecalls_dir)
        self.assertEqual(basecalls.name, "PG1-2")
        self.assertEqual(basecalls.parent, "Rebasecalling")
        self.assertEqual(basecalls.pool, "PG1-2")
        self.assertEqual(basecalls.run, "PG1-4_20240513")
        self.assertEqual(basecalls.bam_pass, str(Path(basecalls_dir).joinpath("pass")))
        self.assertEqual(basecalls.html_report,
                         str(Path(basecalls_dir).joinpath("report_20240513_0829_1A_PAW15419_465bb23f.html")))

class TestHtmlReport(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_html_report(self):
        """
        HtmlReport: load from file
        """
        html_report_file = Path(self.wd).joinpath("report_BLAH.html")
        create_html_report(str(html_report_file))
        html_report = HtmlReport(str(html_report_file))
        self.assertIsNotNone(html_report.extract_json())

class TestBasecallsMetadata(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_basecalls_metadata_no_file(self):
        """
        BasecallsMetadata: no data loaded from file
        """
        data = BasecallsMetadata()
        self.assertEqual(data.flow_cell_id, None)
        self.assertEqual(data.flow_cell_type, None)
        self.assertEqual(data.kit, None)
        self.assertEqual(data.basecalling, None)
        self.assertEqual(data.modifications, None)
        self.assertEqual(data.modified_basecalling, None)
        self.assertEqual(data.trim_barcodes, None)
        self.assertEqual(data.software_versions, None)

    def test_basecalls_metadata_minknow_2_4(self):
        """
        BasecallsMetadata: load from file (MinKNOW 2.4.*)
        """
        html_report_file = Path(self.wd).joinpath("report_BLAH.html")
        create_html_report(str(html_report_file), minknow_version="2.4")
        data = BasecallsMetadata()
        data.load_from_report_html(str(html_report_file))
        self.assertEqual(data.flow_cell_id, "PAW15677")
        self.assertEqual(data.flow_cell_type, "FLO-PRO114M")
        self.assertEqual(data.kit, "SQK-RBK114-24")
        self.assertEqual(data.basecalling, "Super-accurate basecalling, 400 bps")
        self.assertEqual(data.modifications, "5mC & 5hmC")
        self.assertEqual(data.modified_basecalling, "On")
        self.assertEqual(data.trim_barcodes, "Off")
        self.assertEqual(data.software_versions["minknow"], "24.02.19")

    def test_basecalls_metadata_minknow_2_5(self):
        """
        BasecallsMetadata: load from file (MinKNOW 2.5.*)
        """
        html_report_file = Path(self.wd).joinpath("report_BLAH.html")
        create_html_report(str(html_report_file), minknow_version="2.5")
        data = BasecallsMetadata()
        data.load_from_report_html(str(html_report_file))
        self.assertEqual(data.flow_cell_id, "PBC32212")
        self.assertEqual(data.flow_cell_type, "FLO-PRO114M")
        self.assertEqual(data.kit, "SQK-PCB114-24")
        self.assertEqual(data.basecalling, "High-accuracy model 400bps")
        self.assertEqual(data.modifications, None)
        self.assertEqual(data.modified_basecalling, "Off")
        self.assertEqual(data.trim_barcodes, "Off")
        self.assertEqual(data.software_versions["minknow"], "25.03.7")


class TestFlowcellBasecallsInfo(unittest.TestCase):

    def test_flowcell_basecalls_info(self):
        """
        FlowcellBasecallsInfo: load from directory
        """
        pass


class TestIsFlowCellName(unittest.TestCase):

    def test_is_flow_cell_name(self):
        """
        is_flow_cell_name: check flow cell name is identified
        """
        self.assertTrue(is_flow_cell_name("20240422_1039_1E_PAW18123_e180e8f3"))
        self.assertTrue(is_flow_cell_name("20250304_0906_1B_PBC48601_4fb37e4a"))
        self.assertTrue(is_flow_cell_name("20250414_0932_1F_PAY49753_ce3f6499"))
        self.assertFalse(is_flow_cell_name("Rebasecalling"))


class TestGetFlowCellID(unittest.TestCase):

    def test_get_flow_cell_id(self):
        """
        get_flow_cell_id: check ID is extracted
        """
        self.assertEqual(get_flow_cell_id("20240422_1039_1E_PAW18123_e180e8f3"), "PAW18123")
        self.assertEqual(get_flow_cell_id("Rebasecalling"), "")


class TestGetFlowCellDatestamp(unittest.TestCase):

    def test_get_flow_cell_datestamp(self):
        """
        get_flow_cell_datestamp: check datestamp is extracted
        """
        self.assertEqual(get_flow_cell_datestamp("20240422_1039_1E_PAW18123_e180e8f3"), "20240422")
        self.assertEqual(get_flow_cell_datestamp("Rebasecalling"), "")

class TestBarcodeDirs(unittest.TestCase):

    def test_barcode_dirs(self):
        """
        test_barcode_dirs: list directories
        """
        pass
