#!/usr/bin/env python3

# Tests for utils module
import os
import shutil
import tempfile
import unittest
from auto_process_ngs.command import Command
from bcftbx.JobRunner import SimpleJobRunner
from bcf_nanopore.utils import MetadataTabFile
from bcf_nanopore.utils import convert_field_name
from bcf_nanopore.utils import execute_command
from bcf_nanopore.utils import fmt_value
from bcf_nanopore.utils import fmt_yes_no
from pathlib import Path


class TestMetadataTabFile(unittest.TestCase):

    def setUp(self):
        self.metadata_file = tempfile.mkstemp()[1]
        self.projects = list()
        self.lines = list()

    def tearDown(self):
        if self.metadata_file is not None:
            try:
                os.remove(self.metadata_file)
            except PermissionError:
                pass

    def test_empty_metadata_tab_file(self):
        """
        MetaDataTabFile: create and save empty file
        """
        # Make an empty 'file'
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        contents = "#SampleName\tBarcode\tFlowCellID\n"
        self.assertEqual(len(metadata), 0)
        # Save to an actual file and check its contents
        metadata.save(self.metadata_file)
        with open(self.metadata_file, 'rt') as fp:
            self.assertEqual(fp.read(), contents)

    def test_create_new_metadata_tab_file(self):
        """
        MetadataTabFile: create and save MetadataTabFile with content
        """
        # Make new 'file' and add projects
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        self.assertEqual(metadata.index, "SampleName")
        metadata.add_entry('GR1',
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        metadata.add_entry('GR2',
                           barcode="NB04",
                           flow_cell_id="PAW14886")
        contents = "#SampleName\tBarcode\tFlowCellID\nGR1\tNB03\tPAW14886\nGR2\tNB04\tPAW14886\n"
        self.assertEqual(len(metadata), 2)
        # Save to an actual file and check its contents
        metadata.save(self.metadata_file)
        with open(self.metadata_file, 'rt') as fp:
            self.assertEqual(fp.read(), contents)

    def test_read_existing_metadata_tab_file(self):
        """
        MetadataTabFile: read contents from existing MetadataTabFile
        """
        # Create metadata file independently
        data = [dict(SampleName="GR1",
                     Barcode="NB03",
                     FlowCellID="PAW14886"),
                dict(SampleName="GR2",
                     Barcode="NB04",
                     FlowCellID="PAW14886")]
        contents = "#SampleName\tBarcode\tFlowCellID\nGR1\tNB03\tPAW14886\nGR2\tNB04\tPAW14886\n"
        with open(self.metadata_file, 'wt') as fp:
            fp.write(contents)
        # Load and check contents
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"],
                                   filein=self.metadata_file)
        self.assertEqual(len(metadata), 2)
        for x, y in zip(data, metadata):
            for attr in ('SampleName', 'Barcode', 'FlowCellID'):
                self.assertEqual(x[attr], y[attr])

    def test_entries_with_numbers_for_index_values(self):
        """
        MetadataTabFile: handle entries where index values look like numbers
        """
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("1234", barcode="NB03", flow_cell_id="PAW14886")
        metadata.add_entry("56.78", barcode="NB04", flow_cell_id="PAW14886")
        self.assertEqual(metadata[0]['SampleName'], "1234")
        self.assertEqual(metadata[1]['SampleName'], "56.78")

    def test_handle_commented_project(self):
        """
        MetadataTabFile: handle commented entries
        """
        # Create metadata file independently
        data = [dict(SampleName="GR1",
                     Barcode="NB03",
                     FlowCellID="PAW14886"),
                dict(SampleName="#GR2",
                     Barcode="NB04",
                     FlowCellID="PAW14886")]
        contents = "#SampleName\tBarcode\tFlowCellID\nGR1\tNB03\tPAW14886\n#GR2\tNB04\tPAW14886\n"
        with open(self.metadata_file, 'wt') as fp:
            fp.write(contents)
        # Load and check contents
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"],
                                   filein=self.metadata_file)
        self.assertEqual(len(metadata), 2)
        for x, y in zip(data, metadata):
            for attr in ("SampleName", "Barcode", "FlowCellID"):
                self.assertEqual(x[attr], y[attr])

    def test_refuse_to_add_duplicate_index(self):
        """
        MetadataTabFile: refuse to add duplicated index values
        """
        # Make new 'file' and add project
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Attempt to add same project name again
        self.assertRaises(Exception,
                          metadata.add_entry, "GR1",
                          barcode="NB03",
                          flow_cell_id="PAW14886")

    def test_check_entry_is_present(self):
        """
        MetadataTabFile: check that entry is present
        """
        # Make new 'file' and add entry
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Check for existing entry
        self.assertTrue("GR1" in metadata)
        self.assertTrue("#GR1" in metadata)
        # Check for non-existent entry
        self.assertFalse("GR2" in metadata)
        self.assertFalse("#GR2" in metadata)

    def test_check_commented_entry_is_present(self):
        """
        MetadataTabFile: check entry is present when commented
        """
        # Make new 'file' and add entry
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("#GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Check for existing entry
        self.assertTrue("GR1" in metadata)
        self.assertTrue("#GR1" in metadata)
        # Check for non-existent project
        self.assertFalse("GR2" in metadata)
        self.assertFalse("#GR2" in metadata)

    def test_lookup_entry(self):
        """
        MetadataTabFile: look up a specified entry
        """
        # Make new 'file' and add entries
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Fetch data
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")
        # Fetch data for commented name
        data = metadata.lookup_entry("#GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")

    def test_lookup_commented_entry(self):
        """
        MetadataTabFile: check lookup returns specified (commented) entry
        """
        # Make new 'file' and add entries
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("#GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Fetch data
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "#GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")
        # Fetch data for commented name
        data = metadata.lookup_entry("#GR1")
        self.assertEqual(data[0], "#GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")

    def test_lookup_missing_entry_raises_exception(self):
        """
        MetadataTabFile: entry lookup raises KeyError for missing entry
        """
        # Make new 'file'
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        # Check for exceptions
        self.assertRaises(KeyError,
                          metadata.lookup_entry,
                          "GR1")
        self.assertRaises(KeyError,
                          metadata.lookup_entry,
                          "#GR1")

    def test_update_existing_entry(self):
        """
        MetadataTabFile: update data for an existing entry
        """
        # Make new 'file' and add project
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Check initial data
        self.assertTrue("GR1" in metadata)
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")
        # Update flow cell ID value
        metadata.update_entry("GR1", flow_cell_id="PAW14887")
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14887")
        # Update barcode value
        metadata.update_entry("GR1", barcode="NB04")
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB04")
        self.assertEqual(data[2], "PAW14887")

    def test_update_entry_toggles_commenting(self):
        """
        MetadataTabFile: toggle the commenting for an existing entry
        """
        # Make new 'file' and add project
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry("GR1",
                           barcode="NB03",
                           flow_cell_id="PAW14886")
        # Check initial data
        self.assertTrue("GR1" in metadata)
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")
        # Update and comment the entry index
        metadata.update_entry('#GR1')
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")
        # Update and uncomment the entry name
        metadata.update_entry('GR1')
        data = metadata.lookup_entry("GR1")
        self.assertEqual(data[0], "GR1")
        self.assertEqual(data[1], "NB03")
        self.assertEqual(data[2], "PAW14886")

    def test_entries_are_sorted(self):
        """
        MetadataTabFile: entries are sorted into index order
        """
        sorted_contents = """#SampleName\tBarcode\tFlowCellID\nAB1\tNB10\tPAW14886\nXY10\tNB03\tPAW14879\n"""
        metadata = MetadataTabFile(fields=["SampleName",
                                           "Barcode",
                                           "FlowCellID"])
        metadata.add_entry('XY10',
                           barcode="NB03",
                           flow_cell_id="PAW14879")
        metadata.add_entry('AB1',
                           barcode="NB10",
                           flow_cell_id="PAW14886")
        metadata.save(self.metadata_file)
        with open(self.metadata_file, 'rt') as fp:
            self.assertEqual(fp.read(), sorted_contents)


class TestExecuteCommand(unittest.TestCase):

    def setUp(self):
        self.pwd = os.getcwd()
        self.wd = tempfile.mkdtemp()
        os.chdir(self.wd)

    def tearDown(self):
        os.chdir(self.pwd)
        if Path(self.wd).exists():
            shutil.rmtree(self.wd)

    def test_execute_command_using_subprocess(self):
        """
        execute_command: run command using subprocess
        """
        cmd = Command("touch", "test.file")
        self.assertEqual(execute_command(cmd), 0)
        self.assertTrue(Path(self.wd).joinpath("test.file").exists())

    def test_execute_command_using_jobrunner(self):
        """
        execute_command: run command using job runner
        """
        cmd = Command("touch", "test.file")
        runner = SimpleJobRunner(join_logs=True)
        self.assertEqual(execute_command(cmd, runner=runner), 0)
        self.assertTrue(Path(self.wd).joinpath("test.file").exists())


class TestFmtValue(unittest.TestCase):

    def test_fmt_value(self):
        """
        fmt_value: check converted values
        """
        self.assertEqual(fmt_value("hello"), "hello")
        self.assertEqual(fmt_value(""), "")
        self.assertEqual(fmt_value(None), "?")
        self.assertEqual(fmt_value(True), True)
        self.assertEqual(fmt_value(False), False)
        self.assertEqual(fmt_value(0), 0)
        self.assertEqual(fmt_value(1.23), 1.23)


class TestFmtYesNo(unittest.TestCase):

    def test_fmt_yes_no(self):
        """
        fmt_yes_no: check converted values
        """
        self.assertEqual(fmt_yes_no(True), "yes")
        self.assertEqual(fmt_yes_no(False), "yes")
        self.assertEqual(fmt_yes_no(0), "yes")
        self.assertEqual(fmt_yes_no(None), "no")


class TestConvertFieldName(unittest.TestCase):

    def test_convert_field_name(self):
        """
        convert_field_name: check conversions
        """
        self.assertEqual(convert_field_name("sample"), "sample")
        self.assertEqual(convert_field_name("Sample"), "sample")
        self.assertEqual(convert_field_name("SAMPLE"), "sample")
        self.assertEqual(convert_field_name("SampleName"), "sample_name")
        self.assertEqual(convert_field_name("sample_name"), "sample_name")
        self.assertEqual(convert_field_name("Sample Name"), "sample_name")
        self.assertEqual(convert_field_name("sample.name"), "sample_name")
        self.assertEqual(convert_field_name("FlowCellID"), "flow_cell_id")
        self.assertEqual(convert_field_name("flow cell ID"), "flow_cell_id")
