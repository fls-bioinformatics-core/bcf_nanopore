#!/usr/bin/env python3

# Tests for the 'settings' module

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from bcftbx.JobRunner import SimpleJobRunner,GEJobRunner
from bcf_nanopore.settings import Settings


class TestSettings(unittest.TestCase):

    def setUp(self):
        # Create a temp working dir
        self.wd = tempfile.mkdtemp(suffix="TestSettings")

    def tearDown(self):
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_settings_defaults(self):
        """
        Settings: default settings (no .ini file)
        """
        s = Settings()
        # General settings
        self.assertTrue(isinstance(s.general.default_runner,
                                   SimpleJobRunner))
        # QC
        self.assertEqual(s.qc.modules, "nanoplot,fastq_screen")
        self.assertEqual(s.qc.fastq_screen_conf_file, None)
        # Conda
        self.assertFalse(s.conda.enable_conda)
        self.assertEqual(s.conda.env_dir, None)
        # Scheduler
        self.assertEqual(s.scheduler.max_concurrent_jobs, 12)
        self.assertEqual(s.scheduler.max_cores, None)
        self.assertEqual(s.scheduler.poll_interval, 5.0)
        # Runners
        self.assertTrue(isinstance(s.runners.rsync, SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.fastq_screen, SimpleJobRunner))
        self.assertTrue(isinstance(s.runners.nanoplot, SimpleJobRunner))
        # Custom metadata
        self.assertEqual(s.metadata.custom_project_metadata, None)
        self.assertEqual(s.metadata.custom_run_metadata, None)

    def test_settings_from_config_file(self):
        """
        Settings: load settings from config .ini file
        """
        settings_file = Path(self.wd).joinpath("bcf_nanopore.ini")
        with open(settings_file, "wt") as fp:
            fp.write("""[general]
default_runner = SimpleJobRunner

[runners]
rsync = SimpleJobRunner(nslots=4)
""")
        # Load settings
        s = Settings(settings_file)
        # General settings
        self.assertTrue(isinstance(s.general.default_runner,
                                   SimpleJobRunner))
        # Runners
        self.assertTrue(isinstance(s.runners.rsync,
                                   SimpleJobRunner))
        self.assertEqual(s.runners.rsync.nslots, 4)

    def test_settings_custom_metadata_items(self):
        """
        Settings: load custom metadata settings from config .ini file
        """
        settings_file = Path(self.wd).joinpath("bcf_nanopore.ini")
        with open(settings_file, "wt") as fp:
            fp.write("""[metadata]
custom_project_metadata = supplier
custom_run_metadata = order_numbers,analysts
""")
        # Load settings
        s = Settings(settings_file)
        # Custom metadata
        self.assertEqual(s.metadata.custom_project_metadata, "supplier")
        self.assertEqual(s.metadata.custom_run_metadata, "order_numbers,analysts")

    def test_settings_qc_defaults(self):
        """
        Settings: load QC default settings from empty config .ini file
        """
        settings_file = Path(self.wd).joinpath("bcf_nanopore.ini")
        with open(settings_file, "wt") as fp:
            fp.write("""[qc]
""")
        # Load settings
        s = Settings(settings_file)
        # QC settings
        self.assertEqual(s.qc.modules, "nanoplot,fastq_screen")
        self.assertEqual(s.qc.fastq_screen_conf_file, None)

    def test_settings_qc(self):
        """
        Settings: load QC settings from config .ini file
        """
        settings_file = Path(self.wd).joinpath("bcf_nanopore.ini")
        with open(settings_file, "wt") as fp:
            fp.write("""[qc]
modules = nanoplot,fastq_screen,bababa
fastq_screen_conf_file = /path/to/fastqscreen.conf
""")
        # Load settings
        s = Settings(settings_file)
        # QC settings
        self.assertEqual(s.qc.modules, "nanoplot,fastq_screen,bababa")
        self.assertEqual(s.qc.fastq_screen_conf_file, "/path/to/fastqscreen.conf")

    def test_settings_conda(self):
        """
        Settings: load conda settings from config .ini file
        """
        settings_file = Path(self.wd).joinpath("bcf_nanopore.ini")
        with open(settings_file, "wt") as fp:
            fp.write("""[conda]
enable_conda = True
env_dir = /path/to/env_dir
""")
        # Load settings
        s = Settings(settings_file)
        # QC settings
        self.assertTrue(s.conda.enable_conda)
        self.assertEqual(s.conda.env_dir, "/path/to/env_dir")

    def test_settings_reporting_templates(self):
        """
        Settings: load reporting templates from config .ini file
        """
        settings_file = Path(self.wd).joinpath("bcf_nanopore.ini")
        with open(settings_file, "wt") as fp:
            fp.write("""[reporting_templates]
default = name,id,user,pi,application,organism,nsamples,samples
bcf = datestamp,NULL,user,id,#samples,NULL,organism,application,PI,analysis_dir,NULL,primary_data
""")
        # Load settings
        s = Settings(settings_file)
        # Reporting templates
        self.assertEqual([name for name in s.reporting_templates],
                         ["default", "bcf"])
        self.assertEqual(s.reporting_templates.default,
                         "name,id,user,pi,application,organism,"
                         "nsamples,samples")
        self.assertEqual(s.reporting_templates.bcf,
                         "datestamp,NULL,user,id,#samples,NULL,organism,"
                         "application,PI,analysis_dir,NULL,primary_data")
