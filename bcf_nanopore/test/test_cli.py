#!/usr/bin/env python3

# Tests for the 'cli' module

import os
import pwd
import grp
import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.analysis import ProjectAnalysisDir
from bcf_nanopore.mock import MockPromethionDataDir
from bcf_nanopore.mock import MockProjectAnalysisDir
from bcf_nanopore.cli import setup as cli_setup
from bcf_nanopore.cli import fetch as cli_fetch
from bcf_nanopore.cli import report as cli_report

class TestSetupCommand(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_setup_single_run(self):
        """
        setup: create new analysis directory (single run)
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
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
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directory
        run_dir = Path(analysis_dir_path).joinpath("001_PG1-4_20240513")
        self.assertTrue(run_dir.is_dir())
        for f in ["README", "flowcell_basecalls.tsv", "samples.tsv"]:
            self.assertTrue(run_dir.joinpath(f).exists(),
                            f"Expected file {f} not found in run directory")

    def test_setup_multiple_runs(self):
        """
        setup: create new analysis directory (multiple runs)
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-2_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-2_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        data_dir.add_flow_cell("20240529_0830_1A_PAW17328_523ce32d",
                               relpath=Path("PG3-4_20240529").joinpath("PG3-4"))
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
        # Check sub-directories
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        # Check run directories
        for run in ["001_PG1-2_20240513", "002_PG3-4_20240529"]:
            run_dir = Path(analysis_dir_path).joinpath(run)
            self.assertTrue(run_dir.is_dir())
            for f in ["README", "flowcell_basecalls.tsv", "samples.tsv"]:
                self.assertTrue(run_dir.joinpath(f).exists(),
                                f"Expected file {f} not found in run directory "
                                f"'{run}'")

    def test_setup_set_permissions(self):
        """
        setup: create new analysis directory (set permissions)
        """
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        cli_setup(project_dir, "Per Gynt", "Henrik Ibsen", "Methylation study",
                  "Human", top_dir=self.wd, permissions="ugo+rwX")
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
        self.assertFalse(Path(analysis_dir_path).joinpath("samples.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())

    def test_setup_set_group(self):
        """
        setup: create new analysis directory (set group)
        """
        # Find groups for current user
        current_user = pwd.getpwuid(os.getuid()).pw_name
        groups = [g.gr_gid
                  for g in grp.getgrall()
                  if current_user in g.gr_mem]
        if len(groups) < 2:
            raise unittest.SkipTest(f"user '{current_user}' must be in at "
                                    "least two groups for this test")
        # Find a group to set archived files to
        current_gid = os.stat(self.wd).st_gid
        new_group = None
        for gid in groups:
            if gid != current_gid:
                new_group = gid
                break
        self.assertTrue(new_group is not None)
        # Create mock PromethION project
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        project_dir = data_dir.create(self.wd)
        analysis_dir_path = str(Path(self.wd).joinpath("PromethION_Project_001_PerGynt_analysis"))
        cli_setup(project_dir, "Per Gynt", "Henrik Ibsen", "Methylation study",
                  "Human", top_dir=self.wd, group=new_group)
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
        self.assertFalse(Path(analysis_dir_path).joinpath("samples.tsv").exists())
        self.assertTrue(Path(analysis_dir_path).joinpath("ScriptCode").is_dir())
        self.assertTrue(Path(analysis_dir_path).joinpath("logs").is_dir())
        self.assertEqual(Path(analysis_dir_path).stat().st_gid, new_group)


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
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
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
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
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
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        source_dir = os.path.join(self.wd, "source")
        os.mkdir(source_dir)
        project_dir = data_dir.create(source_dir)
        # Fetch subset using job runner
        target_dir = os.path.join(self.wd, "target")
        cli_fetch(project_dir, target_dir, runner="SimpleJobRunner(join_logs=True)")
        self.assertTrue(Path(target_dir).joinpath("PromethION_Project_001_PerGynt").exists())

    def test_fetch_set_permissions(self):
        """
        fetch: copy PromethION data (set permissions)
        """
        # Make source data
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        source_dir = os.path.join(self.wd, "source")
        os.mkdir(source_dir)
        project_dir = data_dir.create(source_dir)
        # Fetch subset
        target_dir = os.path.join(self.wd, "target")
        cli_fetch(project_dir, target_dir, permissions="ugo+rwX")
        self.assertTrue(Path(target_dir).joinpath("PromethION_Project_001_PerGynt").exists())

    def test_fetch_set_group(self):
        """
        fetch: copy PromethION data (set group)
        """
        # Find groups for current user
        current_user = pwd.getpwuid(os.getuid()).pw_name
        groups = [g.gr_gid
                  for g in grp.getgrall()
                  if current_user in g.gr_mem]
        if len(groups) < 2:
            raise unittest.SkipTest(f"user '{current_user}' must be in at "
                                    "least two groups for this test")
        # Find a group to set archived files to
        current_gid = os.stat(self.wd).st_gid
        new_group = None
        for gid in groups:
            if gid != current_gid:
                new_group = gid
                break
        self.assertTrue(new_group is not None)
        # Make source data
        data_dir = MockPromethionDataDir("PromethION_Project_001_PerGynt")
        data_dir.add_flow_cell("20240513_0829_1A_PAW15419_465bb23f",
                               relpath=Path("PG1-4_20240513").joinpath("PG1-2"))
        data_dir.add_basecalls_dir(str(Path("PG1-4_20240513").joinpath("Rebasecalling","PG1-2")),
                                   flow_cell_name="20240513_0829_1A_PAW15419_465bb23f")
        source_dir = os.path.join(self.wd, "source")
        os.mkdir(source_dir)
        project_dir = data_dir.create(source_dir)
        # Fetch subset
        target_dir = os.path.join(self.wd, "target")
        cli_fetch(project_dir, target_dir, group=new_group)
        self.assertTrue(Path(target_dir).joinpath(
            "PromethION_Project_001_PerGynt").exists())
        self.assertEqual(
            Path(target_dir).joinpath("PromethION_Project_001_PerGynt").stat().st_gid, new_group)


class TestReportCommand(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_report_default_summary_mode(self):
        """
        report: default template in 'summary' mode
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
        out_file = Path(self.wd).joinpath("report.txt")
        cli_report(analysis_dir_path, out_file=out_file)
        expected_report = """PromethION_Project_001_PerGynt
==============================
Project name    : PromethION_Project_001_PerGynt
Project ID      : PROMETHION#001


User            : Per Gynt
PI              : Henrik Ibsen
Application     : Methylation study
Organism        : Human

#samples        : 2
samples         : PG1,PG2



"""
        self.assertTrue(out_file.exists())
        with open(out_file, "rt") as fp:
            self.assertEqual(fp.read(), expected_report)

    def test_report_default_tsv_mode(self):
        """
        report: default template in 'tsv' mode
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
        analysis_dir.add_run("PG1-2_20240513",
                             samples={ "PG1": ("NB03", "PAW14589"),
                                       "PG2": ("NB04", "PAW14589")})
        out_file = Path(self.wd).joinpath("report.txt")
        cli_report(analysis_dir_path, mode="tsv", out_file=out_file)
        expected_report = """PromethION_Project_001_PerGynt	PROMETHION#001			Per Gynt	Henrik Ibsen	Methylation study	Human		2	PG1,PG2			
"""
        self.assertTrue(out_file.exists())
        with open(out_file, "rt") as fp:
            self.assertEqual(fp.read(), expected_report)
