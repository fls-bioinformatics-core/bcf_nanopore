#!/usr/bin/env python3

# Tests for the 'analysis' module

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
        # Runners
        self.assertTrue(isinstance(s.runners.rsync,
                                   SimpleJobRunner))
        # Reporting templates

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
