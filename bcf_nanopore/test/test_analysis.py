#!/usr/bin/env python3

# Tests for the 'analysis' module

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.mock import MockPromethionDataDir
from bcf_nanopore.mock import MockProjectAnalysisDir
from bcf_nanopore.analysis import ProjectAnalysisDir
from bcf_nanopore.analysis import FlowcellBasecallsInfo


class TestProjectAnalysisDir(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_project_analysis_dir_create(self):
        """
        ProjectAnalysisDir: create new analysis directory
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human")
        self.assertTrue(analysis_dir.exists())
        self.assertEqual(analysis_dir.path, analysis_dir_path)
        self.assertEqual(analysis_dir.info.name, "PromethION_Project_001_PerGynt")
        self.assertEqual(analysis_dir.info.id, "PROMETHION#001")
        self.assertEqual(analysis_dir.info.platform, "promethion")
        self.assertEqual(analysis_dir.info.data_dir, project_dir)
        self.assertEqual(analysis_dir.info.user, "Per Gynt")
        self.assertEqual(analysis_dir.info.PI, "Henrik Ibsen")
        self.assertEqual(analysis_dir.info.application, "Methylation study")
        self.assertEqual(analysis_dir.info.organism, "Human")
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("basecalling.tsv").exists())
        self.assertFalse(Path(analysis_dir_path).joinpath("samples.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())

    def test_project_analysis_dir_create_with_samples(self):
        """
        ProjectAnalysisDir: create new analysis directory with samples
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human",
                            samples=[("PG1", "NB03", "PAW15419"),
                                     ("PG2", "NB04", "PAW15419"),
                                     ("PG3", "NB05", "PAW15420"),
                                     ("PG4", "NB06", "PAW15420")])
        self.assertTrue(analysis_dir.exists())
        self.assertEqual(analysis_dir.path, analysis_dir_path)
        self.assertEqual(analysis_dir.info.name, "PromethION_Project_001_PerGynt")
        self.assertEqual(analysis_dir.info.id, "PROMETHION#001")
        self.assertEqual(analysis_dir.info.platform, "promethion")
        self.assertEqual(analysis_dir.info.data_dir, project_dir)
        self.assertEqual(analysis_dir.info.user, "Per Gynt")
        self.assertEqual(analysis_dir.info.PI, "Henrik Ibsen")
        self.assertEqual(analysis_dir.info.application, "Methylation study")
        self.assertEqual(analysis_dir.info.organism, "Human")
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("basecalling.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("samples.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())

    def test_project_analysis_dir_load_existing(self):
        """
        ProjectAnalysisDir: load existing analysis project
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir_path = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis").create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            run_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertTrue(analysis_dir.exists())
        self.assertEqual(analysis_dir.path, analysis_dir_path)
        self.assertEqual(analysis_dir.info.name, "PromethION_Project_001_PerGynt")
        self.assertEqual(analysis_dir.info.id, "PROMETHION#001")
        self.assertEqual(analysis_dir.info.platform, "promethion")
        self.assertEqual(analysis_dir.info.data_dir, data_dir)
        self.assertEqual(analysis_dir.info.user, "Per Gynt")
        self.assertEqual(analysis_dir.info.PI, "Henrik Ibsen")
        self.assertEqual(analysis_dir.info.application, "Methylation study")
        self.assertEqual(analysis_dir.info.organism, "Human")
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("basecalling.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("samples.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())

    def test_project_analysis_dir_report(self):
        """
        ProjectAnalysisDir: report analysis project
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir_path = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis").create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            run_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(analysis_dir.report("tsv", "name,#samples,user,pi"),
                         "PromethION_Project_001_PerGynt\t2\tPer Gynt\t"
                         "Henrik Ibsen")


class TestFlowcellBasecallsInfo(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_flowcell_basecalls_info_create_file(self):
        """
        FlowcellBasecallsInfo: create new file
        """
        # New (empty) basecalls info metadata
        basecalls_info = FlowcellBasecallsInfo()
        self.assertEqual(len(basecalls_info), 0)
        # Add data
        basecalls_info.add_base_calls(
            run="Run1",
            pool_name="WT1_WT2_K27CL1_K27CL2",
            sub_dir="WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18",
            flow_cell_id="PBC32212",
            reports="html,json",
            kit="SQK-PCB114-24",
            modifications="none",
            trim_barcodes="Off",
            minknow_version="25.03.7",
            basecalling_model="dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
            file_types="pod5,bam")
        self.assertEqual(len(basecalls_info), 1)
        self.assertEqual(basecalls_info[0]["Run"], "Run1")
        self.assertEqual(basecalls_info[0]["PoolName"], "WT1_WT2_K27CL1_K27CL2")
        self.assertEqual(basecalls_info[0]["SubDir"], "WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18")
        self.assertEqual(basecalls_info[0]["FlowCellID"], "PBC32212")
        self.assertEqual(basecalls_info[0]["Reports"], "html,json")
        self.assertEqual(basecalls_info[0]["Kit"], "SQK-PCB114-24")
        self.assertEqual(basecalls_info[0]["Modifications"], "none")
        self.assertEqual(basecalls_info[0]["TrimBarcodes"], "Off")
        self.assertEqual(basecalls_info[0]["BasecallingModel"], "dna_r10.4.1_e8.2_400bps_hac@v4.3.0")
        self.assertEqual(basecalls_info[0]["FileTypes"], "pod5,bam")
        # Save to file
        basecalls_file = os.path.join(self.wd, "basecalls.tsv")
        self.assertFalse(os.path.exists(basecalls_file))
        basecalls_info.save(basecalls_file)
        self.assertTrue(os.path.exists(basecalls_file))
        # Check contents
        with open(basecalls_file, "rt") as fp:
            contents = fp.read()
            self.assertEqual(contents,
                             """#Run	PoolName	SubDir	FlowCellID	Reports	Kit	Modifications	TrimBarcodes	MinknowVersion	BasecallingModel	FileTypes
Run1	WT1_WT2_K27CL1_K27CL2	WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18	PBC32212	html,json	SQK-PCB114-24	none	Off	25.03.7	dna_r10.4.1_e8.2_400bps_hac@v4.3.0	pod5,bam
""")

    def test_flowcell_basecalls_info_read_from_file(self):
        """
        FlowcellBasecallsInfo: read data from file
        """
        basecalls_file = os.path.join(self.wd, "basecalls.tsv")
        with open(basecalls_file, "wt") as fp:
            fp.write("""#Run	PoolName	SubDir	FlowCellID	Reports	Kit	Modifications	TrimBarcodes	MinknowVersion	BasecallingModel	FileTypes
Run1	WT1_WT2_K27CL1_K27CL2	WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18	PBC32212	html,json	SQK-PCB114-24	none	Off	25.03.7	dna_r10.4.1_e8.2_400bps_hac@v4.3.0	pod5,bam
""")
        basecalls_info = FlowcellBasecallsInfo(basecalls_file)
        self.assertEqual(len(basecalls_info), 1)
        self.assertEqual(basecalls_info[0]["Run"], "Run1")
        self.assertEqual(basecalls_info[0]["PoolName"], "WT1_WT2_K27CL1_K27CL2")
        self.assertEqual(basecalls_info[0]["SubDir"], "WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18")
        self.assertEqual(basecalls_info[0]["FlowCellID"], "PBC32212")
        self.assertEqual(basecalls_info[0]["Reports"], "html,json")
        self.assertEqual(basecalls_info[0]["Kit"], "SQK-PCB114-24")
        self.assertEqual(basecalls_info[0]["Modifications"], "none")
        self.assertEqual(basecalls_info[0]["TrimBarcodes"], "Off")
        self.assertEqual(basecalls_info[0]["BasecallingModel"], "dna_r10.4.1_e8.2_400bps_hac@v4.3.0")
        self.assertEqual(basecalls_info[0]["FileTypes"], "pod5,bam")

    def test_flowcell_basecalls_info_read_from_file_no_minknow_version_or_basecalling_model(self):
        """
        FlowcellBasecallsInfo: read data from file (legacy: no MinKNOW version or basecalling model)
        """
        basecalls_file = os.path.join(self.wd, "basecalls.tsv")
        with open(basecalls_file, "wt") as fp:
            fp.write("""#Run	PoolName	SubDir	FlowCellID	Reports	Kit	Modifications	TrimBarcodes
Run1	WT1_WT2_K27CL1_K27CL2	WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18	PBC32212	yes	SQK-PCB114-24	none	Off
""")
        basecalls_info = FlowcellBasecallsInfo(basecalls_file)
        self.assertEqual(len(basecalls_info), 1)
        self.assertEqual(basecalls_info[0]["Run"], "Run1")
        self.assertEqual(basecalls_info[0]["PoolName"], "WT1_WT2_K27CL1_K27CL2")
        self.assertEqual(basecalls_info[0]["SubDir"], "WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18")
        self.assertEqual(basecalls_info[0]["FlowCellID"], "PBC32212")
        self.assertEqual(basecalls_info[0]["Reports"], "yes")
        self.assertEqual(basecalls_info[0]["Kit"], "SQK-PCB114-24")
        self.assertEqual(basecalls_info[0]["Modifications"], "none")
        self.assertEqual(basecalls_info[0]["TrimBarcodes"], "Off")
        self.assertEqual(basecalls_info[0]["MinknowVersion"], "")
        self.assertEqual(basecalls_info[0]["BasecallingModel"], "")
        self.assertEqual(basecalls_info[0]["FileTypes"], "")
