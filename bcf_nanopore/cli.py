#!/usr/bin/env python3
#
#     cli.py: implement CLI for managing PromethION data for BCF
#     Copyright (C) University of Manchester 2024-2026 Peter Briggs
#

import os
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
from .analysis import RunInfo
from .nanopore.promethion import BasecallsMetadata
from .nanopore.promethion import ProjectDir
from .settings import Settings
from .utils import execute_command
from .utils import fmt_value
from .utils import fmt_yes_no
from . import get_version


# File types
FILE_TYPES = [ "POD5", "FASTQ", "BAM" ]
FILE_TYPES_LOWER = [t.lower() for t in FILE_TYPES]


def get_settings():
    """
    Fetch the configuration settings

    Returns
      Settings: a Settings instance with the default
        configuration parameters loaded
    """
    return Settings()


def get_reporting_templates(settings=None):
    """
    Return the configured reporting templates

    Returns the reporting templates defined in the configuration
    file as a dictionary, with the template names as keys and
    the templates as the associated values.

    Arguments:
      settings (Settings): a Settings instance with the
        configuration parameters to use (otherwise default
        settings will be used)

    Returns:
      Dictionary of reporting templates.
    """
    if settings is None:
        settings = get_settings()
    templates = {}
    for t in [t for t in settings.reporting_templates]:
        templates[t] = settings.reporting_templates[t]
    return templates


def get_custom_metadata_items(settings=None):
    """
    Returns custom project and run metadata items

    Returns the custom metadata items defined in the configuration
    file for project and run levels.

    If no items are defined, returns None; otherwise the items
    are returned as lists.

    Arguments:
      settings (Settings): a Settings instance with the
        configuration parameters to use (otherwise default
        settings will be used)

    Returns:
      Tuple: tuple of (project_metadata, run_metadata) where the
        metadata items are returned as lists, or None if no items
        are defined.
    """
    if settings is None:
        settings = get_settings()
    try:
        project_metadata_items = settings.metadata.custom_project_metadata.split(",")
    except AttributeError:
        project_metadata_items = None
    try:
        run_metadata_items = settings.metadata.custom_run_metadata.split(",")
    except AttributeError:
        run_metadata_items = None
    return (project_metadata_items, run_metadata_items)


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


def extract_metadata(metadata_file, dump_json=False):
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
    custom_project_metadata_items, custom_run_metadata_items = get_custom_metadata_items()
    # Read source project data
    project_name = os.path.basename(os.path.normpath(project_dir))
    # Create analysis dir
    if top_dir is None:
        top_dir = os.getcwd()
    top_dir = os.path.abspath(top_dir)
    analysis_dir = ProjectAnalysisDir(
        os.path.join(top_dir,
                     "%s_analysis" % project_name),
        custom_project_metadata_items=custom_project_metadata_items,
        custom_run_metadata_items=custom_run_metadata_items)
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


def update(path, project_dir, permissions=None, group=None):
    """
    Update an existing analysis directory with new runs

    The analysis directory will be called "<PROJECT>_analysis".

    Arguments:
      path (str): path to PromethION project analysis dir
      project_dir (str): path to PromethION project
      permissions (str): update file permissions on the
        analysis directory using the supplied mode (e.g. 'g+w')
      group (str): update the filesystem group associated
        with the analysis directory to the supplied group name
    """
    custom_project_metadata_items, custom_run_metadata_items = get_custom_metadata_items()
    # Read in data from the analysis directory
    analysis_dir = ProjectAnalysisDir(path,
                                      custom_project_metadata_items=custom_project_metadata_items,
                                      custom_run_metadata_items=custom_run_metadata_items)
    # Do the update
    analysis_dir.update(project_dir)
    # Set permissions and group
    if permissions:
        set_permissions(permissions, analysis_dir.path)
    if group:
        set_group(group, analysis_dir.path)


def metadata(path, items=None):
    """
    Display or update metadata for Promethion project analysis directory

    Arguments:
        path (str): path to PromethION project analysis dir
        items (list): list of metadata update specifications, with
          values of the form "[RUN:]ITEM=VALUE"
    """
    # Get configuration settings
    custom_project_metadata_items, custom_run_metadata_items = get_custom_metadata_items()
    reporting_templates = get_reporting_templates()
    # Read in data
    analysis_dir = ProjectAnalysisDir(path,
                                      custom_project_metadata_items=custom_project_metadata_items,
                                      custom_run_metadata_items=custom_run_metadata_items)
    if items:
        # Update metadata values
        for raw_item in items:
            # Extract run
            try:
                run, item = raw_item.split(":")
            except ValueError:
                run = None
                item = raw_item
            # Split item and value
            try:
                item, value = item.split("=")
            except ValueError:
                raise Exception(f"{raw_item}: invalid metadata item specification")
            # Update the metadata
            if run is None:
                # Update project metadata item
                print(f"...updating '{item}': '{value}'")
                analysis_dir.info[item] = value
                analysis_dir.info.save()
            elif run in analysis_dir.runs:
                # Update run metadata item
                print(f"...updating '{item}' for run '{run}': '{value}'")
                run_info = RunInfo(os.path.join(analysis_dir.path,
                                                analysis_dir.run_dirs[run],
                                                "run.info"),
                                   custom_items=custom_run_metadata_items)
                run_info[item] = value
                run_info.save()
    else:
        # Display project-level metadata
        for item in analysis_dir.info:
            print(f"{item}:\t{analysis_dir.info[item]}")
        # Display run-level metadata for each run
        for run in analysis_dir.runs:
            print(f"\n{run}")
            run_info = RunInfo(os.path.join(analysis_dir.path,
                                            analysis_dir.run_dirs[run],
                                            "run.info"),
                               custom_items=custom_run_metadata_items)
            for item in run_info:
                print(f"\t{item}:\t{run_info[item]}")

def report(path, mode="summary", fields=None, template=None, out_file=None,
           most_recent=None):
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
      most_recent (int): optional, only report this number of most
        recent runs (default is report all runs)
    """
    # Get configuration settings
    custom_project_metadata_items, custom_run_metadata_items = get_custom_metadata_items()
    reporting_templates = get_reporting_templates()
    # Read in data
    analysis_dir = ProjectAnalysisDir(path,
                                      custom_project_metadata_items=custom_project_metadata_items,
                                      custom_run_metadata_items=custom_run_metadata_items)
    # Summary mode
    if mode == "summary":
        report_text = analysis_dir.report_project_summary(most_recent=most_recent)
    elif mode == "runs":
        # Set fields
        if fields is None:
            # Default fields
            fields = "name,id,run,NULL,NULL,user,pi,application,organism,NULL,#samples,samples"
        if template:
            try:
                fields = reporting_templates[template]
            except KeyError:
                raise Exception("%s: undefined template" % template)
        # Report
        report_text = analysis_dir.report_project_runs(fields, most_recent=most_recent)
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
    settings = get_settings()
    default_permissions = settings.general.permissions
    default_group = settings.general.group

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

    # Extract_metadata command
    md_cmd = sp.add_parser("extract_metadata",
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

    # Update command
    update_cmd = sp.add_parser("update",
                              help="Update an analysis directory for a "
                              "Promethion project with new runs")
    update_cmd.add_argument('analysis_dir',
                            help="PromethION analysis directory")
    update_cmd.add_argument('project_dir',
                            help="top level PromethION project directory")
    update_cmd.add_argument('--chmod', action="store",
                            dest="permissions", metavar="PERMISSIONS",
                            default=default_permissions,
                            help="specify permissions for the analysis "
                            "directory. PERMISSIONS should be a string "
                            "recognised by the 'chmod' command (e.g. "
                            "'o-rwX') (default: %s)" %
                            (f"'{default_permissions}'" if default_permissions
                             else "don't set permissions",))
    update_cmd.add_argument('--group', action='store',
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
    report_cmd.add_argument('-r', '--most_recent',
                            metavar="N",
                            action='store',
                            dest='most_recent',
                            default=None,
                            type=int,
                            help="only report N most recent runs")

    # Metadata command
    metadata_cmd = sp.add_parser("metadata",
                                 help="report/update metadata for a PromethION "
                                 "analysis directory")
    metadata_cmd.add_argument('analysis_dir',
                              help="PromethION analysis directory")
    metadata_cmd.add_argument('--set', action="append",
                              dest="item_value", metavar="[RUN:]ITEM=VALUE",
                              help="set metadata ITEM to VALUE; if RUN is "
                              "specified then update the item associated "
                              "with that run, otherwise update the item "
                              "associated with the project")

    # Fetch command
    default_runner = settings.runners.rsync
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
    elif args.command == "update":
        update(args.analysis_dir, args.project_dir,
               permissions=args.permissions, group=args.group)
    elif args.command == "metadata":
        metadata(args.analysis_dir, items=args.item_value)
    elif args.command == "extract_metadata":
        extract_metadata(args.file, dump_json=args.json)
    elif args.command == "report":
        if args.runs:
            mode = "runs"
        else:
            mode = "summary"
        report(args.analysis_dir, mode=mode, fields=args.fields,
               template=args.template, out_file=args.out_file,
               most_recent=args.most_recent)
    elif args.command == "fetch":
        fetch(args.project_dir, args.dest,
              file_types=[x for x in str(args.file_types).split(",")],
              dry_run=args.dry_run, runner=args.runner,
              permissions=args.permissions, group=args.group)
