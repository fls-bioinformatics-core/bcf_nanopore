#!/usr/bin/env python3

# Tests for the 'qc.pipeline' module

import os
import gzip
import shutil
import tempfile
import unittest
from pathlib import Path
from bcf_nanopore.qc.pipeline import NanoporeQC
from bcf_nanopore.mock import MockNanoPlot
from auto_process_ngs.mock import MockSeqtk
from auto_process_ngs.mock import MockFastqScreen


# Polling interval for pipeline
POLL_INTERVAL = 0.1


# Test FASTQ read
FQ_READ = """@432b5bb9-932d-4295-a07f-ab0728acb2e2 runid=50109f28-37e9-4f00-7b64-55ea624ebb57 ch=1234 start_time=2026-11-06T08:23:53.991043+01:00 flow_cell_id=PBC32123 basecall_gpu=NVIDIA_A800_80GB_PCIe protocol_group_id=Example sample_id=Example01 barcode=barcode01 barcode_alias=barcode01 parent_read_id=432b5bb9-932d-4295-a07f-ab0728acb2e2 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.3.0
ACTGTAGTCTACTCGTTGGGGTTGTGCCAGCAATTGGCAAAGTCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTCAGAGGAGTTTAACGTCTCCATGAAATCCCTGTCTGAAGGGCTGACGCATGCTATCATGAGCTCTTGTCGGCCAAAGCAATATCAGCACCAACAGAAACACAAAGACACCGACAACTTTCTTGTGGGCTGCTACCCCAAGGGAATCCTGCAAGAGTCTTGTCCCCAAAGCGCAATCCCAATCCCAATGGTTTGAGCTTAGAGTCTCTGTAAGGGACATGGGTGGCTCTCTT
+
####%%#%$'()($$&''%%%&(&&$$%#$$%&())*)%$#%$%&(+001CCFEKFCDHIKGJSFJQFFEHHHISKSHDGLIJOFGDFA>==AF;8877289<DDE;::;CGA@@?AFIEKHIJKGFELJJIRPGHGFSGFGIPE@??@CCB????DFFH>=11111.))**00:28B??>7>78GB=>>@DA@@ADDHEBLJGS@@HKGSEB,3678:9:=?@@CSIHA:::;<+)()13<;=?SFJSFHLIGKIFLIHG?9421****11,MOSKBLJEEGGIGSHDFFIIKKLHLKOFJADD88:;@><664589AA891--,)*%%
"""


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
            with gzip.open(fq, 'wt') as fp:
                for _ in range(5):
                    fp.write(FQ_READ)
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