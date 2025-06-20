#!/usr/bin/env python3

# Tests for the 'analysis' module

import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.mock import MockPromethionDataDir
from bcf_nanopore.mock import MockProjectAnalysisDir
from bcf_nanopore.analysis import ProjectAnalysisDir


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
        self.assertTrue(Path(analysis_dir_path).joinpath("PromethION_Project_001_PerGynt.tsv").exists())
        self.assertFalse(Path(analysis_dir_path).joinpath("samples.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())

    def test_project_analysis_dir_create_with_samplesheet(self):
        """
        ProjectAnalysisDir: create new analysis directory with sample sheet
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        sample_sheet = Path(self.wd).joinpath("sample_sheet.csv")
        sample_sheet.write_text("""Sample name,Barcode,Flow cell ID
PG1,NB03,PAW15419
PG2,NB04,
PG3,NB05,PAW15420
PG4,NB06,
""")
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human",
                            sample_sheet=str(sample_sheet))
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
        self.assertTrue(Path(analysis_dir_path).joinpath("PromethION_Project_001_PerGynt.tsv").exists())
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
        self.assertTrue(Path(analysis_dir_path).joinpath("PromethION_Project_001_PerGynt.tsv").exists())
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
