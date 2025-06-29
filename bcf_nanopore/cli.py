#!/usr/bin/env python3
#
#     cli.py: implement CLI for managing PromethION data for BCF
#     Copyright (C) University of Manchester 2024-2025 Peter Briggs
#

import os
import json
from argparse import ArgumentParser
from auto_process_ngs.command import Command
from bcftbx.JobRunner import fetch_runner
from .analysis import ProjectAnalysisDir
from .nanopore.promethion import BasecallsMetadata
from .nanopore.promethion import ProjectDir
from .utils import execute_command
from .utils import fmt_value
from .utils import fmt_yes_no


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
                     "TrimBarcodes"]))
    for fc in project.flow_cells:
        kit = fmt_value(fc.metadata.kit)
        modifications = ("none" if fc.metadata.modified_basecalling == "Off"
                         else fmt_value(fc.metadata.modifications))
        trim_barcodes = fmt_value(fc.metadata.trim_barcodes)
        has_report = fmt_yes_no(fc.html_report)
        print('\t'.join([str(s)for s in (fc.run,
                                         fc.pool,
                                         fc,
                                         fc.id,
                                         has_report,
                                         kit,
                                         modifications,
                                         trim_barcodes)]))
    for bc in project.basecalls_dirs:
        flow_cell_id = fmt_value(bc.metadata.flow_cell_id)
        kit = fmt_value(bc.metadata.kit)
        modifications = ("none" if bc.metadata.modified_basecalling == "Off"
                         else fmt_value(bc.metadata.modifications))
        trim_barcodes = fmt_value(bc.metadata.trim_barcodes)
        has_report = fmt_yes_no(bc.html_report)
        print('\t'.join([str(s)for s in (bc.run,
                                         bc.name,
                                         bc,
                                         flow_cell_id,
                                         has_report,
                                         kit,
                                         modifications,
                                         trim_barcodes)]))


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
    if metadata_file.endswith(".html"):
        data = BasecallsMetadata()
        data.load_from_report_html(metadata_file)
        if dump_json:
            print(data.json())
        else:
            print("Flow cell ID         : %s" % data.flow_cell_id)
            print("Flow cell type       : %s" % data.flow_cell_type)
            print("Kit type             : %s" % data.kit)
            print("Modified basecalling : %s" % data.modified_basecalling)
            print("Modified base context: %s" % data.modifications)
            print("Barcode trimming     : %s" % data.trim_barcodes)
    elif metadata_file.endswith(".json"):
        with open(metadata_file, "rt") as fp:
            data = json.load(fp)
        try:
            print(json.dumps(data, sort_keys=True, indent=4))
        except BrokenPipeError:
            pass


def setup(project_dir, user, PI, application=None, organism=None,
          samples_csv=None, top_dir=None):
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


def report(path, mode="summary", fields=None, template=None, out_file=None):
    """
    Report on Promethion project analysis directory
    """
    # Templates
    _templates = {
        # Default: for spreadsheet
        'default':
        "name,id,NULL,NULL,user,pi,application,organism,NULL,"
        "nsamples,samples,NULL,NULL,NULL",
        # BCF: for downstream spreadsheet
        'bcf': "datestamp,NULL,user,id,#samples,NULL,organism,application,PI,analysis_dir,NULL,primary_data",
        # Summary: for reporting run for downstream analysis
        'summary': "name,id,datestamp,platform,analysis_dir,NULL,"
        "user,pi,application,organism,primary_data,comments",
    }
    # Read in data
    analysis_dir = ProjectAnalysisDir(path)
    # Set fields
    if fields is None:
        if template is None:
            if mode == "summary":
                template = "summary"
            else:
                template = "default"
        try:
            fields = _templates[template]
        except KeyError:
            raise Exception("%s: undefined template" % template)
    # Report
    report_text = analysis_dir.report(mode, fields)
    if out_file is None:
        print(report_text)
    else:
        print(f"Writing to {out_file}")
        with open(out_file, 'wt') as fp:
            fp.write(report_text + '\n')


def fetch(project_dir, target_dir, dry_run=False, runner=None):
    """
    Fetch the BAM files and reports for a Promethion run

    Arguments:
      project_dir (str): path to PromethION project
        (can local or remote)
      target_dir (str): path to top-level directory to
        copy project files into
      dry_run (bool): if True then do dry run rsync only
        (default is to actually fetch the data)
      runner (str): job runner definition to use to
        execute the fetch operations
    """
    # Clean the project dir path
    project_dir = project_dir.rstrip(os.sep)
    # Project name
    project_name = os.path.basename(project_dir)
    # Fetch job runner
    if runner is not None:
        runner = fetch_runner(runner)
        print(f"Using job runner '{runner}'")
    # Example rsync to only fetch BAM and index files:
    # rsync --dry-run -av -m --include="*/" \
    # --include="bam_pass/*/*.bam" --include="bam_pass/*/*.bai" \
    # --include="pass/*/*.bam" --include="pass/*/*.bai" \
    # --exclude="*" \
    # <PromethION_PROJECT_DIR> .
    rsync_bams = Command('rsync')
    if dry_run:
        rsync_bams.add_args('--dry-run')
    rsync_bams.add_args('-av',
                        '-m',
                        '--include=*/',
                        '--include=bam_pass/*/*.bam',
                        '--include=bam_pass/*/*.bai',
                        '--include=pass/*/*.bam',
                        '--include=pass/*/*.bai',
                        '--exclude=*',
                        project_dir,
                        target_dir)
    print("Transferring BAM files with command: %s" % rsync_bams)
    status = execute_command(rsync_bams, runner=runner)
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
                           '--exclude=*',
                           project_dir,
                           target_dir)
    print("Transferring report files with command: %s" % rsync_reports)
    status = execute_command(rsync_reports, runner=runner)
    if status != 0:
        raise Exception("fetch: failed to transfer reports")


def bcf_nanopore_main():

    # Main parser
    p = ArgumentParser()
    sp = p.add_subparsers(dest='command')

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
                        "items (HTML reports only)")

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
                            choices=['default', 'bcf', 'summary'],
                            help="specify template used to set fields "
                            "for reporting")
    report_cmd.add_argument('-f', '--fields',
                            help="specify fields to report (comma-separated "
                            "list; overrides template from '-t')")
    report_cmd.add_argument('-o', '--out_file',
                            help="write summary to specified file rather "
                            "than stdout")

    # Fetch command
    fetch_cmd = sp.add_parser("fetch",
                              help="fetch BAM files from PromethION project "
                              "directory")
    fetch_cmd.add_argument('project_dir',
                           help="top level PromethION project directory")
    fetch_cmd.add_argument('dest',
                           help="destination directory (copy of top-level "
                           "directory will be created under this)")
    fetch_cmd.add_argument('--dry-run', action="store_true",
                           help="dry run only (no data will be copied)")
    fetch_cmd.add_argument('-r', '--runner', action="store",
                           help="job runner to use (optional)")

    # Process command line
    args = p.parse_args()

    # Execute command
    if args.command == "info":
        info(args.project_dir)
    elif args.command == "setup":
        setup(args.project_dir, user=args.user, PI=args.pi,
              application=args.application, organism=args.organism,
              samples_csv=args.samples_csv, top_dir=args.parent_dir)
    elif args.command == "metadata":
        metadata(args.file, dump_json=args.json)
    elif args.command == "report":
        report(args.analysis_dir, mode=args.mode, fields=args.fields,
               template=args.template, out_file=args.out_file)
    elif args.command == "fetch":
        fetch(args.project_dir, args.dest, dry_run=args.dry_run,
              runner=args.runner)
