#!/usr/bin/env python3
#
#     cli.py: implement CLI for managing PromethION data for BCF
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs
#

import os
import json
import tempfile
import shutil
from argparse import ArgumentParser
from auto_process_ngs.command import Command
from auto_process_ngs.fileops import copy
from auto_process_ngs.fileops import set_group
from auto_process_ngs.fileops import set_permissions
from bcftbx.JobRunner import fetch_runner
from bcftbx.JobRunner import BaseJobRunner
from .analysis import ProjectAnalysisDir
from .nanopore.promethion import BasecallsMetadata
from .nanopore.promethion import ProjectDir
from .settings import Settings
from .utils import execute_command
from .utils import fmt_value
from .utils import fmt_yes_no

# Configuration
__settings = Settings()

# File types
FILE_TYPES = [ "POD5", "FASTQ", "BAM" ]
FILE_TYPES_LOWER = [t.lower() for t in FILE_TYPES]

# Reporting templates
REPORTING_TEMPLATES = {
    # Default: for spreadsheet
    'default': "name,id,NULL,NULL,user,pi,application,organism,NULL,"
    "nsamples,samples,NULL,NULL,NULL",
    # BCF: for downstream spreadsheet
    'bcf': "datestamp,NULL,user,id,#samples,NULL,organism,application,PI,analysis_dir,NULL,primary_data",
    # Summary: for reporting run for downstream analysis
    'summary': "name,id,datestamp,platform,analysis_dir,NULL,"
    "user,pi,application,organism,primary_data,comments",
}
for t in [t for t in __settings.reporting_templates]:
    REPORTING_TEMPLATES[t] = __settings.reporting_templates[t]

def config():
    """
    Print configuration information
    """
    settings = Settings(resolve_undefined=True)
    print(settings.report_settings(exclude_undefined=False))


def info(project_dir):
    """
    Print information on a Promethion project directory

    Outputs TSV data with descriptions of runs, pools etc.
    under a project.

    Arguments:
      project_dir (str): path to the top-level PromethION
        project directory
    """
    # Top level
    project = ProjectDir(os.path.abspath(project_dir))
    # Make a metadata file
    print('\t'.join(["#Run",
                     "PoolName",
                     "SubDir",
                     "FlowCellID",
                     "Reports",
                     "Kit",
                     "Modifications",
                     "TrimBarcodes",
                     "MinKNOWVersion",
                     "BasecallingModel",
                     "FileTypes"]))
    for fc in project.flow_cells:
        run = ("-" if fc.run is None else fc.run)
        kit = fmt_value(fc.metadata.kit)
        modifications = ("none" if fc.metadata.modified_basecalling == "Off"
                         else fmt_value(fc.metadata.modifications))
        trim_barcodes = fmt_value(fc.metadata.trim_barcodes)
        try:
            minknow_version = fc.metadata.software_versions["minknow"]
        except (TypeError, KeyError):
            minknow_version = "?"
        reports = fc.report_types
        if reports:
            reports = ",".join(reports)
        else:
            reports = "none"
        file_types = fc.file_types
        if file_types:
            file_types = ",".join(file_types)
        else:
            file_types = "none"
        basecalling_model = fc.metadata.basecalling_model
        if basecalling_model is None:
            basecalling_model = fc.metadata.basecalling_config
        basecalling_model = fmt_value(basecalling_model)
        print('\t'.join([str(s) for s in (run,
                                          fc.pool,
                                          os.path.relpath(fc.path,
                                                          project.path),
                                          fc.id,
                                          reports,
                                          kit,
                                          modifications,
                                          trim_barcodes,
                                          minknow_version,
                                          basecalling_model,
                                          file_types)]))
    for bc in project.basecalls_dirs:
        flow_cell_id = fmt_value(bc.metadata.flow_cell_id)
        run = ("-" if bc.run is None else bc.run)
        if bc.pool:
            pool = bc.pool
        else:
            pool = bc.name
        kit = fmt_value(bc.metadata.kit)
        modifications = ("none" if bc.metadata.modified_basecalling == "Off"
                         else fmt_value(bc.metadata.modifications))
        trim_barcodes = fmt_value(bc.metadata.trim_barcodes)
        try:
            minknow_version = bc.metadata.software_versions["minknow"]
        except (TypeError, KeyError):
            minknow_version = "?"
        reports = bc.report_types
        if reports:
            reports = ",".join(reports)
        else:
            reports = "none"
        file_types = bc.file_types
        if file_types:
            file_types = ",".join(file_types)
        else:
            file_types = "none"
        basecalling_model = bc.metadata.basecalling_model
        if basecalling_model is None:
            basecalling_model = bc.metadata.basecalling_config
        basecalling_model = fmt_value(basecalling_model)
        print('\t'.join([str(s) for s in (run,
                                          pool,
                                          os.path.relpath(bc.path,
                                                          project.path),
                                          flow_cell_id,
                                          reports,
                                          kit,
                                          modifications,
                                          trim_barcodes,
                                          minknow_version,
                                          basecalling_model,
                                          file_types)]))


def metadata(metadata_file, dump_json=False):
    """
    Extract metadata from HTML report or JSON file

    Arguments:
      metadata_file (str): path to HTML report or JSON
        file to extract data from
      dump_json (bool): for HTML reports, if False
        (default) then only print extract metadata items;
        otherwise dump the extract JSON data
    """
    data = BasecallsMetadata()
    data.load_from_report(metadata_file)
    if metadata_file.endswith(".html"):
        if dump_json:
            try:
                print(data.html_json())
            except BrokenPipeError:
                pass
        else:
            print("Flow cell ID         : %s" % data.flow_cell_id)
            print("Flow cell type       : %s" % data.flow_cell_type)
            print("Kit type             : %s" % data.kit)
            print("Modified basecalling : %s" % data.modified_basecalling)
            print("Modified base context: %s" % data.modifications)
            print("Barcode trimming     : %s" % data.trim_barcodes)
            print("Software versions    : %s" % data.software_versions)
    elif metadata_file.endswith(".json"):
        if dump_json:
            try:
                print(data.json())
            except BrokenPipeError:
                pass
        else:
            print("Basecalling model    : %s" % data.basecalling_model)
            print("Basecalling config   : %s" % data.basecalling_config)


def setup(project_dir, user, PI, application=None, organism=None,
          samples_csv=None, top_dir=None, permissions=None, group=None):
    """
    Set up a new analysis directory for a Promethion project

    The analysis directory will be called "<PROJECT>_analysis".

    Information about the samples can be supplied via a CSV
    file with the format:

    ::

        Header line
        SAMPLE,BARCODE[,FLOWCELL]

    Arguments:
      project_dir (str): path to PromethION project
      user (str): user(s) associated with the project
      PI (str): principal investigator(s)
      application (str): experimental application(s)
      organism (str): associated origanism(s)
      samples_csv (str): path to CSV file with sample information
      top_dir (str): directory to make analysis directory
        under (defaults to current directory)
      permissions (str): update file permissions on the
        analysis directory using the supplied mode (e.g. 'g+w')
      group (str): update the filesystem group associated
        with the analysis directory to the supplied group name
    """
    # Read source project data
    project_name = os.path.basename(os.path.normpath(project_dir))
    # Fetch sample information
    samples = []
    if samples_csv:
        # Samples data should be a CSV format file with an
        # initial header line (which is ignored) followed by
        # lines with either 2 or 3 fields
        prev_flowcell = None
        ignore_line = True
        with open(samples_csv, "rt") as fp:
            for line in fp:
                if ignore_line:
                    # Ignore first line
                    ignore_line = False
                    continue
                try:
                    # Three fields: sample, barcode, flowcell ID
                    sample, barcode, flowcell = line.strip().split(',')
                except ValueError:
                    # Two fields: sample, barcode
                    # Assumes same flow cell as previous line
                    sample, barcode = line.strip().split(',')
                if flowcell:
                    prev_flowcell = flowcell
                else:
                    flowcell = prev_flowcell
                samples.append((sample, barcode, flowcell))
    # Create analysis dir
    if top_dir is None:
        top_dir = os.getcwd()
    top_dir = os.path.abspath(top_dir)
    analysis_dir = ProjectAnalysisDir(
        os.path.join(top_dir,
                     "%s_analysis" % project_name))
    analysis_dir.create(project_dir,
                        user=user,
                        PI=PI,
                        application=application,
                        organism=organism,
                        samples=samples)
    # Set permissions and group
    if permissions:
        set_permissions(permissions, analysis_dir.path)
    if group:
        set_group(group, analysis_dir.path)


def report(path, mode="summary", fields=None, template=None, out_file=None):
    """
    Report on Promethion project analysis directory

    Arguments:
      path (str): path to PromethION project analysis dir
      mode (str): either "summary" or "tsv" (default: summary)
      fields (str): optional, list of fields to report (overrides
        "template")
      template (str): optional, name of a pre-defined template
        (i.e. set of fields) to use from configuration
      out_file (str): optional, file to write the report to
        (default is to write to stdout)
    """
    # Read in data
    analysis_dir = ProjectAnalysisDir(path)
    # Set fields
    if fields is None:
        if template is None:
            template = "default"
        try:
            fields = REPORTING_TEMPLATES[template]
        except KeyError:
            raise Exception("%s: undefined template" % template)
    # Report
    report_text = analysis_dir.report(mode, fields)
    if out_file:
        # File extension for report file
        ext = "tsv" if mode == "tsv" else "txt"
        # Temporary copy
        temp_dir = tempfile.mkdtemp()
        temp_file = os.path.join(temp_dir,
                                 f"{analysis_dir.info.id}.{ext}")
        with open(temp_file, "wt") as fp:
            fp.write(report_text + '\n')
        # Copy to final location
        copy(temp_file, str(out_file))
        shutil.rmtree(temp_dir)
        print(f"Written to {out_file}")
    else:
        print(report_text)


def fetch(project_dir, target_dir, file_types=None, dry_run=False,
          runner=None, permissions=None, group=None):
    """
    Fetch the BAM files and reports for a Promethion run

    Arguments:
      project_dir (str): path to PromethION project
        (can local or remote)
      target_dir (str): path to top-level directory to
        copy project files into
      file_types (list): list of file types to include
        (can be one or more of "pod5", "fastq", "bam";
        defaults to "bam" if none supplied)
      dry_run (bool): if True then do dry run rsync only
        (default is to actually fetch the data)
      runner (str): job runner definition to use to
        execute the fetch operations
      permissions (str): update file permissions on the
        copied files and directories using the supplied
        mode (e.g. 'g+w')
      group (str): update the filesystem group associated
        with the copied files to the supplied group name
    """
    # Clean the project dir path
    project_dir = project_dir.rstrip(os.sep)
    # Project name
    project_name = os.path.basename(project_dir)
    # Fetch job runner
    if runner is not None and not isinstance(runner, BaseJobRunner):
            runner = fetch_runner(runner)
    print(f"Using job runner '{runner}'")
    # File types to include
    if not file_types:
        file_types = ["bam"]
    for t in file_types:
        if t.lower() not in FILE_TYPES_LOWER:
            raise Exception(f"fetch: unknown file type requested: "
                            f"'{t}'")
    file_types = [t.lower() for t in file_types]
    # Example rsync to only fetch BAM and index files:
    # rsync --dry-run -av -m --include="*/" \
    # --include="bam_pass/*/*.bam" --include="bam_pass/*/*.bai" \
    # --include="pass/*/*.bam" --include="pass/*/*.bai" \
    # --exclude="*" \
    # <PromethION_PROJECT_DIR> .
    rsync_data = Command('rsync')
    if dry_run:
        rsync_data.add_args("--dry-run")
    rsync_data.add_args("-av",
                        "-m",
                        "--include=*/")
    if "pod5" in file_types:
        rsync_data.add_args(
            "--include=pod5/*/*.pod5")
    if "fastq" in file_types:
        rsync_data.add_args(
            "--include=fastq_pass/*/*.fastq",
            "--include=fastq_pass/*/*.fastq.gz",
            "--include=pass/*/*.fastq",
            "--include=pass/*/*.fastq.gz")
    if "bam" in file_types:
        rsync_data.add_args(
            "--include=bam_pass/*/*.bam",
            "--include=bam_pass/*/*.bai",
            "--include=pass/*/*.bam",
            "--include=pass/*/*.bai")
    rsync_data.add_args("--exclude=*")
    if permissions:
        rsync_data.add_args(f"--chmod={permissions}")
    rsync_data.add_args(project_dir,
                        target_dir)
    print("Transferring data files with command: %s" % rsync_data)
    status = execute_command(rsync_data, runner=runner)
    if status != 0:
        raise Exception("fetch: failed to transfer data")
    # Example to fetch reports and sample sheets:
    # rsync --dry-run -av -m --include="*/" --include="report_*" \
    # --include="sample_sheet_*" --exclude="*" \
    # <PromethION_PROJECT_DIR> .
    rsync_reports = Command('rsync')
    if dry_run:
        rsync_reports.add_args('--dry-run')
    rsync_reports.add_args('-av',
                           '-m',
                           '--include=*/',
                           '--include=report_*',
                           '--include=sample_sheet_*',
                           '--exclude=*')
    if permissions:
        rsync_reports.add_args(f"--chmod={permissions}")
    rsync_reports.add_args(project_dir,
                           target_dir)
    print("Transferring report files with command: %s" % rsync_reports)
    status = execute_command(rsync_reports, runner=runner)
    if status != 0:
        raise Exception("fetch: failed to transfer reports")
    # Set the group
    if group is not None:
        print(f"Setting group on copied files to '{group}'")
        if not dry_run:
            set_group(group, os.path.join(target_dir, project_name))


def bcf_nanopore_main():

    # Defaults
    default_permissions = __settings.general.permissions
    default_group = __settings.general.group

    # Main parser
    p = ArgumentParser()
    sp = p.add_subparsers(dest='command')

    # Config command
    config_cmd = sp.add_parser("config",
                               help="Print configuration information")

    # Info command
    info_cmd = sp.add_parser("info",
                             help="Get information on a Promethion project "
                             "directory")
    info_cmd.add_argument('project_dir',
                          help="top level PromethION project directory")

    # Metadata command
    md_cmd = sp.add_parser("metadata",
                           help="Extract metadata from HTML or JSON report")
    md_cmd.add_argument('file',
                        help="HTML or JSON report file")
    md_cmd.add_argument('-j', '--json', action='store_true',
                        help="dump extracted JSON data instead of metadata "
                        "items")

    # Setup command
    setup_cmd = sp.add_parser("setup",
                              help="Set up a new analysis directory for a "
                              "Promethion project")
    setup_cmd.add_argument('project_dir',
                           help="top level PromethION project directory")
    setup_cmd.add_argument('parent_dir', nargs="?",
                           help="create project analysis directory under "
                           "PARENT_DIR (defaults to current directory)")
    setup_cmd.add_argument('-u', '--user', action='store', required=True,
                           help="User associated with project")
    setup_cmd.add_argument('-p', '--pi', action='store', required=True,
                           help="PI associated with project")
    setup_cmd.add_argument('-a', '--application', action='store',
                           help="applications(s) associated with project")
    setup_cmd.add_argument('-o', '--organism', action='store',
                           help="organism(s) associated with project")
    setup_cmd.add_argument('-s', '--samples_csv', action='store',
                           help="CSV file with 'sample,barcode[,flowcell]' "
                           "information")
    setup_cmd.add_argument('--chmod', action="store",
                           dest="permissions", metavar="PERMISSIONS",
                           default=default_permissions,
                           help="specify permissions for the analysis "
                           "directory. PERMISSIONS should be a string "
                           "recognised by the 'chmod' command (e.g. "
                           "'o-rwX') (default: %s)" %
                           (f"'{default_permissions}'" if default_permissions
                            else "don't set permissions",))
    setup_cmd.add_argument('--group', action='store',
                           default=default_group,
                           help="specify the name of group for the "
                           "analysis directory (default: %s)" %
                           (f"'{default_group}'" if default_group
                            else "don't set group",))

    # Report command
    report_cmd = sp.add_parser("report",
                               help="report metadata from a PromethION "
                               "analysis directory")
    report_cmd.add_argument('analysis_dir',
                            help="PromethION analysis directory")
    report_cmd.add_argument('-m', '--mode',
                            choices=['summary', 'tsv'], default='summary',
                            help="specify reporting mode")
    report_cmd.add_argument('-t', '--template',
                            choices=[t for t in REPORTING_TEMPLATES],
                            help="specify template used to set fields "
                            "for reporting")
    report_cmd.add_argument('-f', '--fields',
                            help="specify fields to report (comma-separated "
                            "list; overrides template from '-t')")
    report_cmd.add_argument('-o', '--out_file',
                            help="write summary to specified file rather "
                            "than stdout")

    # Fetch command
    default_runner = __settings.runners.rsync
    fetch_cmd = sp.add_parser("fetch",
                              help="fetch BAM files from PromethION project "
                              "directory")
    fetch_cmd.add_argument('project_dir',
                           help="top level PromethION project directory")
    fetch_cmd.add_argument('dest',
                           help="destination directory (copy of top-level "
                           "directory will be created under this)")
    fetch_cmd.add_argument('--files', action="store",
                           dest="file_types", metavar="FILETYPES",
                           default="bam",
                           help="specify types of data files to copy "
                           "as a comma-separated list (e.g. 'fastq,bams') "
                           "(valid types are 'pod5', 'fastq', 'bams'; "
                           "default: 'bam')")
    fetch_cmd.add_argument('--chmod', action="store",
                           dest="permissions", metavar="PERMISSIONS",
                           default=default_permissions,
                           help="specify permissions for the copied files. "
                           "PERMISSIONS should be a string recognised by the "
                           "'chmod' command (e.g. 'o-rwX') (default: %s)" %
                           (f"'{default_permissions}'" if default_permissions
                            else "don't set permissions",))
    fetch_cmd.add_argument('--group', action='store',
                           default=default_group,
                           help="specify the name of group for the copied "
                           "files (default: %s)" %
                           (f"'{default_group}'" if default_group
                            else "don't set group",))
    fetch_cmd.add_argument('--dry-run', action="store_true",
                           help="dry run only (no data will be copied)")
    fetch_cmd.add_argument('-r', '--runner', action="store",
                           default=default_runner,
                           help=f"job runner to use (default: "
                           f"'{default_runner}')")

    # Process command line
    args = p.parse_args()

    # Execute command
    if args.command == "config":
        config()
    if args.command == "info":
        info(args.project_dir)
    elif args.command == "setup":
        setup(args.project_dir, user=args.user, PI=args.pi,
              application=args.application, organism=args.organism,
              samples_csv=args.samples_csv, top_dir=args.parent_dir,
              permissions=args.permissions, group=args.group)
    elif args.command == "metadata":
        metadata(args.file, dump_json=args.json)
    elif args.command == "report":
        report(args.analysis_dir, mode=args.mode, fields=args.fields,
               template=args.template, out_file=args.out_file)
    elif args.command == "fetch":
        fetch(args.project_dir, args.dest,
              file_types=[x for x in str(args.file_types).split(",")],
              dry_run=args.dry_run, runner=args.runner,
              permissions=args.permissions, group=args.group)
