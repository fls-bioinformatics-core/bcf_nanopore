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
from . import get_version

# Configuration
__settings = Settings()

# File types
FILE_TYPES = [ "POD5", "FASTQ", "BAM" ]
FILE_TYPES_LOWER = [t.lower() for t in FILE_TYPES]

# Reporting templates from configuration
REPORTING_TEMPLATES = {}
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
                     "SubDir",
                     "FlowCellID",
                     "Reports",
                     "Kit",
                     "Modifications",
                     "TrimBarcodes",
                     "MinKNOWVersion",
                     "BasecallingModel",
                     "FileTypes"]))
    for run in project.runs:
        for fc in run.flow_cells:
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
            print('\t'.join([str(s) for s in (run.name,
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
        for bc in run.basecalls_dirs:
            flow_cell_id = fmt_value(bc.metadata.flow_cell_id)
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
            print('\t'.join([str(s) for s in (run.name,
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


def setup(project_dir, user, PI, application=None, organism=None, top_dir=None,
          permissions=None, group=None):
    """
    Set up a new analysis directory for a Promethion project

    The analysis directory will be called "<PROJECT>_analysis".

    Arguments:
      project_dir (str): path to PromethION project
      user (str): user(s) associated with the project
      PI (str): principal investigator(s)
      application (str): experimental application(s)
      organism (str): associated origanism(s)
      top_dir (str): directory to make analysis directory
        under (defaults to current directory)
      permissions (str): update file permissions on the
        analysis directory using the supplied mode (e.g. 'g+w')
      group (str): update the filesystem group associated
        with the analysis directory to the supplied group name
    """
    # Read source project data
    project_name = os.path.basename(os.path.normpath(project_dir))
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
                        organism=organism)
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
      mode (str): either "summary" or "runs" (default: summary)
      fields (str): optional, list of fields to report in "runs"
        mode (when it overrides "template", if supplied);
        otherwise it is ignored
      template (str): optional, name of a pre-defined template
        (i.e. set of fields) to use from configuration in
        "runs" mode (will be overridden by "fields", if
        supplied); otherwise it is ignored
      out_file (str): optional, file to write the report to
        (default is to write to stdout)
    """
    # Read in data
    analysis_dir = ProjectAnalysisDir(path)
    # Summary mode
    if mode == "summary":
        report_text = analysis_dir.report_project_summary()
    elif mode == "runs":
        # Set fields
        if fields is None:
            # Default fields
            fields = "name,id,run,NULL,NULL,user,pi,application,organism,NULL,#samples,samples"
        if template:
            try:
                fields = REPORTING_TEMPLATES[template]
            except KeyError:
                raise Exception("%s: undefined template" % template)
        # Report
        report_text = analysis_dir.report_project_runs(fields)
    else:
        raise Exception("%s: unknown reporting mode" % mode)
    if out_file:
        # File extension for report file
        ext = "tsv" if mode == "runs" else "txt"
        # Temporary copy
        temp_dir = tempfile.mkdtemp()
        temp_file = os.path.join(temp_dir,
                                 f"{analysis_dir.info.id.lower()}.{ext}")
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
    # --include="bam_pass/**.bam" --include="bam_pass/**.bai" \
    # --include="pass/**.bam" --include="pass/**.bai" \
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
            "--include=pod5/**.pod5")
    if "fastq" in file_types:
        rsync_data.add_args(
            "--include=fastq_pass/**.fastq",
            "--include=fastq_pass/**.fastq.gz",
            "--include=pass/**.fastq",
            "--include=pass/**.fastq.gz")
    if "bam" in file_types:
        rsync_data.add_args(
            "--include=bam_pass/**.bam",
            "--include=bam_pass/**.bai",
            "--include=pass/**.bam",
            "--include=pass/**.bai")
    rsync_data.add_args("--exclude=*")
    if permissions:
        rsync_data.add_args(f"--chmod={permissions}")
    rsync_data.add_args(project_dir,
                        target_dir)
    print("Transferring data files with command: %s" % rsync_data)
    status = execute_command(rsync_data, runner=runner)
    if status != 0:
        raise Exception("fetch: failed to transfer data")
    # Example to fetch reports, sample sheets and summaries:
    # rsync --dry-run -av -m --include="*/" --include="report_*" \
    # --include="sample_sheet_*" --include="sequencing_summary*.txt" \
    # --exclude="*" <PromethION_PROJECT_DIR> .
    rsync_reports = Command('rsync')
    if dry_run:
        rsync_reports.add_args('--dry-run')
    rsync_reports.add_args('-av',
                           '-m',
                           '--include=*/',
                           '--include=report_*',
                           '--include=sample_sheet_*',
                           '--include=sequencing_summary*.txt',
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

    # Version
    p.add_argument('--version', action='version',
                   version=("%%(prog)s %s" % get_version()))

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
    mutex = report_cmd.add_mutually_exclusive_group()
    mutex.add_argument('--summary',
                       action='store_true',
                       dest='summary',
                       default=False,
                       help="print summary report suitable for informaticians "
                       "(default mode)")
    mutex.add_argument('--runs',action='store_true',dest='runs',
                       default=False,
                       help="print tab-delimited line (one per run) "
                       "suitable for injection into a spreadsheet")
    report_cmd.add_argument('-f', '--fields',
                            action='store',
                            dest='fields',
                            default=None,
                            help="fields to report in --runs mode (overrides "
                            "fields specified by --template)")
    report_cmd.add_argument('-t', '--template',
                            action='store',
                            dest='template',
                            default=None,
                            help="name of template with fields to report in "
                            "--runs mode (templates should be defined in the "
                            "config file)")
    report_cmd.add_argument('--file',
                            action='store',
                            dest='out_file',
                            default=None,
                            help="write report to OUT_FILE; destination can be "
                            "a local file, or a remote file specified as "
                            "[[USER@]HOST:]PATH (default is to write to stdout)")

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
                           "as a comma-separated list (e.g. 'fastq,bam') "
                           "(valid types are 'pod5', 'fastq', 'bam'; "
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
              top_dir=args.parent_dir, permissions=args.permissions,
              group=args.group)
    elif args.command == "metadata":
        metadata(args.file, dump_json=args.json)
    elif args.command == "report":
        if args.runs:
            mode = "runs"
        else:
            mode = "summary"
        report(args.analysis_dir, mode=mode, fields=args.fields,
               template=args.template, out_file=args.out_file)
    elif args.command == "fetch":
        fetch(args.project_dir, args.dest,
              file_types=[x for x in str(args.file_types).split(",")],
              dry_run=args.dry_run, runner=args.runner,
              permissions=args.permissions, group=args.group)
