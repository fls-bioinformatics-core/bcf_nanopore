#!/usr/bin/env python3

# Tests for the 'cli' module

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.analysis import ProjectAnalysisDir
from bcf_nanopore.mock import MockPromethionDataDir
from bcf_nanopore.cli import setup as cli_setup
from bcf_nanopore.cli import fetch as cli_fetch

class TestSetupCommand(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_setup(self):
        """
        setup: create new analysis directory
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        cli_setup(project_dir, "Per Gynt", "Henrik Ibsen", "Methylation study",
                  "Human", top_dir=self.wd)
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
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

    def test_setup_with_samplesheet(self):
        """
        setup: create new analysis directory with samples CSV file
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        samples_csv = Path(self.wd).joinpath("samples.csv")
        samples_csv.write_text("""Sample name,Barcode,Flow cell ID
PG1,NB03,PAW15419
PG2,NB04,
PG3,NB05,PAW15420
PG4,NB06,
""")
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        cli_setup(project_dir, "Per Gynt", "Henrik Ibsen", "Methylation study",
                  "Human", samples_csv=str(samples_csv), top_dir=self.wd)
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
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


class TestFetchCommand(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_fetch(self):
        """
        fetch: copy PromethION data
        """
        # Make source data
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        source_dir = os.path.join(self.wd, "source")
        os.mkdir(source_dir)
        project_dir = data_dir.create(source_dir)
        # Fetch subset
        target_dir = os.path.join(self.wd, "target")
        cli_fetch(project_dir, target_dir)
        self.assertTrue(Path(target_dir).joinpath("PromethION_Project_001_PerGynt").exists())

    def test_fetch_trailing_slash_on_source_name(self):
        """
        fetch: copy PromethION data (trailing slash on source name)
        """
        # Make source data
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        source_dir = os.path.join(self.wd, "source")
        os.mkdir(source_dir)
        project_dir = data_dir.create(source_dir)
        # Fetch subset
        target_dir = os.path.join(self.wd, "target")
        cli_fetch(project_dir + os.sep, target_dir)
        self.assertTrue(Path(target_dir).joinpath("PromethION_Project_001_PerGynt").exists())

    def test_fetch_using_jobrunner(self):
        """
        fetch: copy PromethION data using job runner
        """
        # Make source data
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f", run="PG1-4_20240513", pool="PG1-2")
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        source_dir = os.path.join(self.wd, "source")
        os.mkdir(source_dir)
        project_dir = data_dir.create(source_dir)
        # Fetch subset using job runner
        target_dir = os.path.join(self.wd, "target")
        cli_fetch(project_dir, target_dir, runner="SimpleJobRunner(join_logs=True)")
        self.assertTrue(Path(target_dir).joinpath("PromethION_Project_001_PerGynt").exists())
