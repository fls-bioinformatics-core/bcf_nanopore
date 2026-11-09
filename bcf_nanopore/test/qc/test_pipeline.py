#!/usr/bin/env python3

# Tests for the 'qc.pipeline' module

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.qc.pipeline import NanoporeQC
from bcf_nanopore.mock import MockNanoPlot
from bcf_nanopore.mock import mock_fastq
from auto_process_ngs.mock import MockSeqtk
from auto_process_ngs.mock import MockFastqScreen


# Polling interval for pipeline
POLL_INTERVAL = 0.1


class TestNanoporeQC(unittest.TestCase):

    def setUp(self):
        # Temporary working directory for the test
        self.wd = tempfile.mkdtemp()
        # Directory for mock executables
        self.bin = os.path.join(self.wd,"bin")
        os.mkdir(self.bin)
        # Update PATH to include mock executables
        self.path = os.environ["PATH"]
        os.environ['PATH'] = "%s:%s" % (self.bin, self.path)
        # Data directory within the temporary working directory
        self.data = os.path.join(self.wd, "data")
        os.mkdir(self.data)
        # Add empty FastqScreen conf file
        self.fastq_screen_conf_file = os.path.join(
            self.data, "fastq_screen.conf")
        with open(self.fastq_screen_conf_file, 'wt') as fp:
            fp.write("")

    def tearDown(self):
        # Restore PATH
        os.environ["PATH"] = self.path
        # Remove working directory
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_nanopore_qc(self):
        """
        qc.pipeline: test the Nanopore QC pipeline (default modules)
        """
        # Make mock QC executables
        MockSeqtk.create(os.path.join(self.bin, "seqtk"))
        MockFastqScreen.create(os.path.join(self.bin, "fastq_screen"))
        MockNanoPlot.create(os.path.join(self.bin, "NanoPlot"))
        # Create mock Fastqs
        fqs = ["PBC321213_pass_barcode01_45678a10_5432b345_0.fastq.gz",
               "PBC321213_pass_barcode01_45678a10_5432b345_1.fastq.gz"]
        fastqs = []
        for fq in fqs:
            fq = os.path.join(self.data, fq)
            fastqs.append(fq)
            mock_fastq(fq, nreads=5)
        # Set up and run the QC
        nanopore_qc = NanoporeQC()
        out_dir = os.path.join(self.wd, "nanopore_qc")
        status = nanopore_qc.run(
            fastqs,
            self.fastq_screen_conf_file,
            out_dir,
            filter_fastqs_min_reads=1,
            poll_interval=POLL_INTERVAL,
            max_jobs=1)
        # Check the outputs
        self.assertEqual(status,0)