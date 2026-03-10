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

    def test_project_analysis_dir_create_single_run(self):
        """
        ProjectAnalysisDir: create new analysis directory (single run)
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
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
        # Check top-level analysis directory
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
        self.assertEqual(analysis_dir.info.runs, "PG1-4_20240513")
        self.assertEqual(analysis_dir.runs, ["PG1-4_20240513"])
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directory
        run_dir = Path(analysis_dir_path).joinpath("001_PG1-4_20240513")
        self.assertTrue(run_dir.is_dir())
        for f in ["README", "flowcell_basecalls.tsv", "samples.tsv", "run.info"]:
            self.assertTrue(run_dir.joinpath(f).exists(),
                            f"Expected file {f} not found in run directory")
        # Check datestamps
        self.assertEqual(analysis_dir.datestamp(), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG1-4_20240513"), "20240513")
        self.assertEqual(analysis_dir.datestamp_short(), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG1-4_20240513"), "240513")

    def test_project_analysis_dir_create_multiple_runs(self):
        """
        ProjectAnalysisDir: create new analysis directory (multiple runs)
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-2_20240513").joinpath("PG1-2"))
        data_dir.add_flow_cell("20240529_0830_1A_PAW17328_523ce32d",
                               relpath=Path("PG3-4_20240529").joinpath("PG3-4"))
        data_dir.add_basecalls_dir(str(Path("PG3-4_20240529").joinpath("Rebasecalling","PG3-4")),
                                   flow_cell_name="20240529_0830_1A_PAW17328_523ce32d")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human")
        # Check top-level analysis directory
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
        self.assertEqual(analysis_dir.info.runs, "PG1-2_20240513,PG3-4_20240529")
        self.assertEqual(analysis_dir.runs, ["PG1-2_20240513", "PG3-4_20240529"])
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directories
        for run in ["001_PG1-2_20240513", "002_PG3-4_20240529"]:
            run_dir = Path(analysis_dir_path).joinpath(run)
            self.assertTrue(run_dir.is_dir())
            for f in ["README", "flowcell_basecalls.tsv", "samples.tsv", "run.info"]:
                self.assertTrue(run_dir.joinpath(f).exists(),
                                f"Expected file {f} not found in run directory "
                                f"'{run}'")
        # Check datestamps
        self.assertEqual(analysis_dir.datestamp(), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG1-2_20240513"), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG3-4_20240529"), "20240529")
        self.assertEqual(analysis_dir.datestamp_short(), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG1-2_20240513"), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG3-4_20240529"), "240529")

    def test_project_analysis_dir_create_custom_project_metadata(self):
        """
        ProjectAnalysisDir: create new analysis directory with custom project metadata
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path,
                                          custom_project_metadata_items=["order_numbers", "analyst"])
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human")
        # Check top-level analysis directory
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
        self.assertEqual(analysis_dir.info.runs, "PG1-4_20240513")
        self.assertEqual(analysis_dir.info.order_numbers, None)
        self.assertEqual(analysis_dir.info.analyst, None)
        self.assertEqual(analysis_dir.runs, ["PG1-4_20240513"])
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directory
        run_dir = Path(analysis_dir_path).joinpath("001_PG1-4_20240513")
        self.assertTrue(run_dir.is_dir())
        for f in ["README", "flowcell_basecalls.tsv", "samples.tsv", "run.info"]:
            self.assertTrue(run_dir.joinpath(f).exists(),
                            f"Expected file {f} not found in run directory")
        # Check run metadata items
        expected_items = set(["Run name"])
        with open(run_dir.joinpath("run.info"), "rt") as fp:
            for line in fp:
                item = line.split("\t")[0]
                self.assertTrue(item in expected_items, f"{item}: unexpected item")
                expected_items.remove(item)
        self.assertEqual(len(expected_items), 0,
                         f"{expected_items}: items not found in 'run.info'")
        # Check datestamps
        self.assertEqual(analysis_dir.datestamp(), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG1-4_20240513"), "20240513")
        self.assertEqual(analysis_dir.datestamp_short(), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG1-4_20240513"), "240513")

    def test_project_analysis_dir_create_custom_run_metadata(self):
        """
        ProjectAnalysisDir: create new analysis directory with custom run metadata
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path,
                                          custom_run_metadata_items=["order_numbers", "analyst"])
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human")
        # Check top-level analysis directory
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
        self.assertEqual(analysis_dir.info.runs, "PG1-4_20240513")
        self.assertEqual(analysis_dir.runs, ["PG1-4_20240513"])
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directory
        run_dir = Path(analysis_dir_path).joinpath("001_PG1-4_20240513")
        self.assertTrue(run_dir.is_dir())
        for f in ["README", "flowcell_basecalls.tsv", "samples.tsv", "run.info"]:
            self.assertTrue(run_dir.joinpath(f).exists(),
                            f"Expected file {f} not found in run directory")
        # Check run metadata items
        expected_items = set(["Run name", "Order numbers", "Analyst"])
        with open(run_dir.joinpath("run.info"), "rt") as fp:
            for line in fp:
                item = line.split("\t")[0]
                self.assertTrue(item in expected_items, f"{item}: unexpected item")
                expected_items.remove(item)
        self.assertEqual(len(expected_items), 0,
                         f"{expected_items}: items not found in 'run.info'")
        # Check datestamps
        self.assertEqual(analysis_dir.datestamp(), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG1-4_20240513"), "20240513")
        self.assertEqual(analysis_dir.datestamp_short(), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG1-4_20240513"), "240513")

    def test_project_analysis_dir_update_with_new_runs(self):
        """
        ProjectAnalysisDir: update analysis directory with new runs
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-2_20240513").joinpath("PG1-2"))
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertFalse(analysis_dir.exists())
        analysis_dir.create(project_dir,
                            user="Per Gynt",
                            PI="Henrik Ibsen",
                            application="Methylation study",
                            organism="Human")
        # Check top-level analysis directory
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
        self.assertEqual(analysis_dir.info.runs, "PG1-2_20240513")
        self.assertEqual(analysis_dir.runs, ["PG1-2_20240513"])
        self.assertTrue(Path(analysis_dir_path).joinpath("README").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("project.info").exists())
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directories
        for run in ["001_PG1-2_20240513"]:
            run_dir = Path(analysis_dir_path).joinpath(run)
            self.assertTrue(run_dir.is_dir())
            for f in ["README", "flowcell_basecalls.tsv", "samples.tsv", "run.info"]:
                self.assertTrue(run_dir.joinpath(f).exists(),
                                f"Expected file {f} not found in run directory "
                                f"'{run}'")
        # Check datestamps
        self.assertEqual(analysis_dir.datestamp(), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG1-2_20240513"), "20240513")
        self.assertEqual(analysis_dir.datestamp_short(), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG1-2_20240513"), "240513")
        # Add new run with extra flowcell and basecalls
        data_dir.add_flow_cell("20240529_0830_1A_PAW17328_523ce32d",
                               relpath=Path("PG3-4_20240529").joinpath("PG3-4"))
        data_dir.add_basecalls_dir(str(Path("PG3-4_20240529").joinpath("Rebasecalling","PG3-4")),
                                   flow_cell_name="20240529_0830_1A_PAW17328_523ce32d")
        data_dir.update(self.wd)
        # Update the analysis directory
        analysis_dir.update(project_dir)
        # Check runs metadata
        self.assertEqual(analysis_dir.info.runs, "PG1-2_20240513,PG3-4_20240529")
        self.assertEqual(analysis_dir.runs, ["PG1-2_20240513", "PG3-4_20240529"])
        # Check run directories
        for run in ["001_PG1-2_20240513", "002_PG3-4_20240529"]:
            run_dir = Path(analysis_dir_path).joinpath(run)
            self.assertTrue(run_dir.is_dir())
            for f in ["README", "flowcell_basecalls.tsv", "samples.tsv", "run.info"]:
                self.assertTrue(run_dir.joinpath(f).exists(),
                                f"Expected file {f} not found in run directory "
                                f"'{run}'")
        # Check datestamps
        self.assertEqual(analysis_dir.datestamp(), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG1-2_20240513"), "20240513")
        self.assertEqual(analysis_dir.datestamp("PG3-4_20240529"), "20240529")
        self.assertEqual(analysis_dir.datestamp_short(), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG1-2_20240513"), "240513")
        self.assertEqual(analysis_dir.datestamp_short("PG3-4_20240529"), "240529")

    def test_project_analysis_dir_load_existing(self):
        """
        ProjectAnalysisDir: load existing analysis project
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
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
        self.assertEqual(analysis_dir.info.runs, "PG1-2_20240513,PG3-4_20240529")
        self.assertEqual(analysis_dir.runs, ["PG1-2_20240513", "PG3-4_20240529"])

    def test_project_analysis_dir_load_existing_custom_project_metadata(self):
        """
        ProjectAnalysisDir: load existing analysis project with custom project metadata
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001",
            extra_project_metadata={ "Order numbers": "#00123,#00456",
                                     "Analyst": "Ben Franklin" })
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
        self.assertEqual(analysis_dir.info.runs, "PG1-2_20240513,PG3-4_20240529")
        self.assertEqual(analysis_dir.info.order_numbers, "#00123,#00456")
        self.assertEqual(analysis_dir.info.analyst, "Ben Franklin")
        self.assertEqual(analysis_dir.runs, ["PG1-2_20240513", "PG3-4_20240529"])

    def test_project_analysis_dir_load_existing_legacy_run_dir_names(self):
        """
        ProjectAnalysisDir: load existing analysis project (legacy run dir names)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        # Rename run dirs to legacy names
        os.rename(os.path.join(analysis_dir_path, "001_PG1-2_20240513"),
                  os.path.join(analysis_dir_path, "01-PG1-2_20240513"))
        os.rename(os.path.join(analysis_dir_path, "002_PG3-4_20240529"),
                  os.path.join(analysis_dir_path, "02-PG3-4_20240529"))
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
        self.assertEqual(analysis_dir.info.runs, "PG1-2_20240513,PG3-4_20240529")
        self.assertEqual(analysis_dir.runs, ["PG1-2_20240513", "PG3-4_20240529"])

    def test_project_analysis_dir_single_run_report_project_summary(self):
        """
        ProjectAnalysisDir: report summary for analysis project (single run)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(analysis_dir.report_project_summary(),
                         "PromethION_Project_001_PerGynt\n"
                         "==============================\n"
                         "\n"
                         "Project name    : PromethION_Project_001_PerGynt\n"
                         "Project ID      : PROMETHION#001\n"
                         "User            : Per Gynt\n"
                         "PI              : Henrik Ibsen\n"
                         "Application     : Methylation study\n"
                         "Organism        : Human\n"
                         f"Analysis dir    : {analysis_dir_path}\n"
                         "\n"
                         "This project has 1 run:\n"
                         "\n"
                         "- PG1-2_20240513:\t2 samples (PG1, PG2)")

    def test_project_analysis_dir_single_run_report_project_runs(self):
        """
        ProjectAnalysisDir: report runs for analysis project (single run)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={"PG1": ("NB03", "PAW14589"),
                                      "PG2": ("NB04", "PAW14589")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs("name,run,#samples,user,pi"),
            "PromethION_Project_001_PerGynt\tPG1-2_20240513\t2\tPer Gynt\tHenrik Ibsen")

    def test_project_analysis_dir_multiple_runs_report_project_summary(self):
        """
        ProjectAnalysisDir: report summary for analysis project (multiple runs)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(analysis_dir.report_project_summary(),
                         "PromethION_Project_001_PerGynt\n"
                         "==============================\n"
                         "\n"
                         "Project name    : PromethION_Project_001_PerGynt\n"
                         "Project ID      : PROMETHION#001\n"
                         "User            : Per Gynt\n"
                         "PI              : Henrik Ibsen\n"
                         "Application     : Methylation study\n"
                         "Organism        : Human\n"
                         f"Analysis dir    : {analysis_dir_path}\n"
                         "\n"
                         "This project has 2 runs:\n"
                         "\n"
                         "- PG1-2_20240513:	2 samples (PG1, PG2)\n"
                         "- PG3-4_20240529:	2 samples (PG3, PG4)")

    def test_project_analysis_dir_multiple_runs_report_project_summary_most_recent(self):
        """
        ProjectAnalysisDir: report summary for analysis project (multiple runs)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir.add_run("PG5-6_20240602",
                             samples={ "PG5": ("NB11", "PAW16072"),
                                       "PG6": ("NB12", "PAW16072")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(analysis_dir.report_project_summary(most_recent=1),
                         "PromethION_Project_001_PerGynt\n"
                         "==============================\n"
                         "\n"
                         "Project name    : PromethION_Project_001_PerGynt\n"
                         "Project ID      : PROMETHION#001\n"
                         "User            : Per Gynt\n"
                         "PI              : Henrik Ibsen\n"
                         "Application     : Methylation study\n"
                         "Organism        : Human\n"
                         f"Analysis dir    : {analysis_dir_path}\n"
                         "\n"
                         "This project has 1 new run:\n"
                         "\n"
                         "- PG5-6_20240602:	2 samples (PG5, PG6)\n"
                         "\n"
                         "This is in addition to 2 previous runs")

    def test_project_analysis_dir_multiple_runs_report_project_runs(self):
        """
        ProjectAnalysisDir: report runs for analysis project (multiple runs)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs("name,run,#samples,user,pi"),
            "PromethION_Project_001_PerGynt\tPG1-2_20240513\t2\tPer Gynt\tHenrik Ibsen\n"
            "PromethION_Project_001_PerGynt\tPG3-4_20240529\t2\tPer Gynt\tHenrik Ibsen")

    def test_project_analysis_dir_multiple_runs_report_project_runs_most_recent(self):
        """
        ProjectAnalysisDir: report most recent runs for analysis project
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir.add_run("PG5-6_20240602",
                             samples={ "PG5": ("NB11", "PAW16072"),
                                       "PG6": ("NB12", "PAW16072")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs("name,run,#samples,user,pi", most_recent=1),
            "PromethION_Project_001_PerGynt\tPG5-6_20240602\t2\tPer Gynt\tHenrik Ibsen")

    def test_project_analysis_dir_no_samples_report_project_summary(self):
        """
        ProjectAnalysisDir: report summary for analysis project (no samples)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513")
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(analysis_dir.report_project_summary(),
                         "PromethION_Project_001_PerGynt\n"
                         "==============================\n"
                         "\n"
                         "Project name    : PromethION_Project_001_PerGynt\n"
                         "Project ID      : PROMETHION#001\n"
                         "User            : Per Gynt\n"
                         "PI              : Henrik Ibsen\n"
                         "Application     : Methylation study\n"
                         "Organism        : Human\n"
                         f"Analysis dir    : {analysis_dir_path}\n"
                         "\n"
                         "This project has 1 run:\n"
                         "\n"
                         "- PG1-2_20240513:\tno samples")

    def test_project_analysis_dir_no_samples_run_report_project_runs(self):
        """
        ProjectAnalysisDir: report runs for analysis project (no samples)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513")
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs("name,run,#samples,user,pi"),
            "PromethION_Project_001_PerGynt\tPG1-2_20240513\t0\tPer Gynt\tHenrik Ibsen")

    def test_project_analysis_dir_report_project_runs_composite_fields(self):
        """
        ProjectAnalysisDir: report runs for analysis project (composite fields)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs("[_]:name+run,#samples+sample_names,[/]:pi+user"),
            "PromethION_Project_001_PerGynt_PG1-2_20240513\t2 PG1,PG2\tHenrik Ibsen/Per Gynt\n"
            "PromethION_Project_001_PerGynt_PG3-4_20240529\t2 PG3,PG4\tHenrik Ibsen/Per Gynt")

    def test_project_analysis_dir_report_project_runs_bad_composite_fields(self):
        """
        ProjectAnalysisDir: report runs for analysis project (composite fields)
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis")
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        analysis_dir.add_run("PG3-4_20240529",
                             samples={ "PG3": ("NB07", "PAW15894"),
                                       "PG4": ("NB08", "PAW15894")})
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertRaises(ValueError,
                          analysis_dir.report_project_runs,
                          "[_:name+run,#samples")
        self.assertRaises(ValueError,
                          analysis_dir.report_project_runs,
                          "[_]name+run,#samples")
        self.assertRaises(KeyError,
                          analysis_dir.report_project_runs,
                          "_]:name+run,#samples")

    def test_project_analysis_dir_single_run_report_project_runs_custom_metadata(self):
        """
        ProjectAnalysisDir: report runs for analysis project with custom metadata
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis",)
        analysis_dir.add_run("PG1-2_20240513",
                             samples={"PG1": ("NB03", "PAW14589"),
                                      "PG2": ("NB04", "PAW14589")},
                             metadata={ "Order numbers": "#00123" })
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001",
            extra_project_metadata={ "Analysts": "Sam Beckett" })
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs("name,run,#samples,user,pi,analysts,order_numbers"),
            "PromethION_Project_001_PerGynt\tPG1-2_20240513\t2\tPer Gynt\tHenrik Ibsen\tSam Beckett\t#00123")

    def test_project_analysis_dir_single_run_report_project_runs_datestamps(self):
        """
        ProjectAnalysisDir: handle datestamps when reporting runs
        """
        data_dir = "/mnt/data/PromethION_Project_001_PerGynt"
        analysis_dir = MockProjectAnalysisDir("PromethION_Project_001_PerGynt_analysis",)
        analysis_dir.add_run("PG1-2_20240513",
                             samples={"PG1": ("NB03", "PAW14589"),
                                      "PG2": ("NB04", "PAW14589")},)
        analysis_dir_path = analysis_dir.create(
            self.wd,
            user="Per Gynt",
            principal_investigator="Henrik Ibsen",
            application="Methylation study",
            organism="Human",
            data_dir=data_dir,
            project_id="PROMETHION#001")
        analysis_dir = ProjectAnalysisDir(analysis_dir_path)
        self.assertEqual(
            analysis_dir.report_project_runs(
                "datestamp,datestamp_short,run_datestamp,run_datestamp_short,name,run,#samples"),
            "20240513\t240513\t20240513\t240513\tPromethION_Project_001_PerGynt\tPG1-2_20240513\t2")


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
                             """#Run	SubDir	FlowCellID	Reports	Kit	Modifications	TrimBarcodes	MinknowVersion	BasecallingModel	FileTypes
Run1	WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18	PBC32212	html,json	SQK-PCB114-24	none	Off	25.03.7	dna_r10.4.1_e8.2_400bps_hac@v4.3.0	pod5,bam
""")

    def test_flowcell_basecalls_info_read_from_file(self):
        """
        FlowcellBasecallsInfo: read data from file
        """
        basecalls_file = os.path.join(self.wd, "basecalls.tsv")
        with open(basecalls_file, "wt") as fp:
            fp.write("""#Run	SubDir	FlowCellID	Reports	Kit	Modifications	TrimBarcodes	MinknowVersion	BasecallingModel	FileTypes
Run1	WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18	PBC32212	html,json	SQK-PCB114-24	none	Off	25.03.7	dna_r10.4.1_e8.2_400bps_hac@v4.3.0	pod5,bam
""")
        basecalls_info = FlowcellBasecallsInfo(basecalls_file)
        self.assertEqual(len(basecalls_info), 1)
        self.assertEqual(basecalls_info[0]["Run"], "Run1")
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
            fp.write("""#Run	SubDir	FlowCellID	Reports	Kit	Modifications	TrimBarcodes
Run1	WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18	PBC32212	yes	SQK-PCB114-24	none	Off
""")
        basecalls_info = FlowcellBasecallsInfo(basecalls_file)
        self.assertEqual(len(basecalls_info), 1)
        self.assertEqual(basecalls_info[0]["Run"], "Run1")
        self.assertEqual(basecalls_info[0]["SubDir"], "WT1_WT2_K27CL1_K27CL2/20250616_0817_1F_PBC32212_40107e18")
        self.assertEqual(basecalls_info[0]["FlowCellID"], "PBC32212")
        self.assertEqual(basecalls_info[0]["Reports"], "yes")
        self.assertEqual(basecalls_info[0]["Kit"], "SQK-PCB114-24")
        self.assertEqual(basecalls_info[0]["Modifications"], "none")
        self.assertEqual(basecalls_info[0]["TrimBarcodes"], "Off")
        self.assertEqual(basecalls_info[0]["MinknowVersion"], "")
        self.assertEqual(basecalls_info[0]["BasecallingModel"], "")
        self.assertEqual(basecalls_info[0]["FileTypes"], "")
