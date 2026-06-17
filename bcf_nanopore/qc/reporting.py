#!/usr/bin/env python3

# Reporting functions for Nanopore QC


import os
import json
import textwrap
from .utils import split_fastq_name
from ..nanopore.promethion import ProjectDir
from auto_process_ngs.docwriter import Document
from auto_process_ngs.docwriter import Table
from auto_process_ngs.docwriter import Img
from auto_process_ngs.docwriter import Link
from auto_process_ngs.qc.plots import uscreenplot
from auto_process_ngs.utils import ZipMaker


def report_nanopore_qc(out_dir, project_dir=None):
    """
    Generate HTML report
    """
    # Collect the QC outputs
    qc_outputs = {}

    # Read counts
    read_counts_dir = os.path.join(out_dir, "read_counts")
    if os.path.exists(read_counts_dir):
        counts_files = sorted([os.path.join(read_counts_dir, f)
                               for f in os.listdir(read_counts_dir)
                               if f.endswith("_count.out")])
        if counts_files:
            qc_outputs["read_counts"] = {}
            for f in counts_files:
                name = os.path.basename(f)[:-len("_count.out")]
                with open(f, "rt") as fp:
                    try:
                        nreads = int(fp.read().strip())
                    except Exception as ex:
                        print(f"{f}: failed to extract read count: {ex}")
                    qc_outputs["read_counts"][name] = nreads

    # NanoPlot outputs
    nanoplot_dir = os.path.join(out_dir, "nanoplot")
    if os.path.exists(nanoplot_dir):
        qc_outputs["nanoplot"] = {}
        fq_groups = os.listdir(nanoplot_dir)
        for name in fq_groups:
            nanoplot_grp_dir = os.path.join(nanoplot_dir, name)
            nanoplot_report = os.path.join(nanoplot_grp_dir,
                                           "NanoPlot-report.html")
            if os.path.exists(nanoplot_report):
                qc_outputs["nanoplot"][name] = nanoplot_grp_dir

    # FastqScreen output
    screens_dir = os.path.join(out_dir, "fastq_screen")
    if os.path.exists(screens_dir):
        qc_outputs["fastq_screen"] = {}
        screen_files = sorted([f for f in os.listdir(screens_dir)
                               if f.endswith(".txt") and "_screen_" in f])
        try:
            # Attempt to resort using trailing index
            screen_files = sorted(
                screen_files,
                key=lambda f:
                int(f.split(".")[0].split("_screen_")[0].split("_")[-1]))
        except ValueError:
            pass
        for f in screen_files:
            fq_base = f.split(".")[0].split("_screen_")[0]
            fc, bc, _ = split_fastq_name(f)
            name = fc
            if bc:
                name += "_" + bc
            qc_outputs["fastq_screen"][name] = f

    # Collect flowcell/barcode combinations
    fc_bc_names = set()
    for metric in ("read_counts",
                   "nanoplot",
                   "fastq_screen"):
        if metric in qc_outputs:
            fc_bc_names.update(
                set([x for x in qc_outputs[metric].keys()]))
    fc_bc_names = sorted(list(fc_bc_names))

    # Collect project and run information
    if project_dir:
        project = ProjectDir(project_dir)
        project_name = project.name
        runs = {}
        for run in project.runs:
            run_name = run.name
            runs[run_name] = []
            for flowcell in run.flow_cells:
                runs[run_name].append(flowcell.id)
    else:
        project_name = None
        runs = {}
        flowcells = set()
        for name in fc_bc_names:
            try:
                fc, _ = name.split("_")
            except ValueError:
                fc = name
            flowcells.add(fc)
        runs[None] = sorted(list(flowcells))

    # Dump raw output data as JSON
    with open(os.path.join(out_dir, "nanopore_qc.json"), "wt") as fp:
        json.dump({"project_dir": project_dir,
                   "fc_bc_names": fc_bc_names,
                   "runs": runs,
                   "qc_outputs": qc_outputs},
                  fp, sort_keys=True, indent=4)

    # Initialise report
    title = "NanoporeQC report"
    if project_name:
        title += f": {project_name}"
    report = Document(title)
    report.add_css_rule(textwrap.dedent("""
        html { font-family: sans-serif; }
        h1 {
          background-color: #42AEC2;
          color: white;
        }
        h2 {
          background-color: #8CC63F;
          color: white;
          display: inline-block;
        }
        table.summary {
          border: solid 1px grey;
          background-color: white;
          margin: 10 10;
          font-size: 80%;
        }
        table.summary th {
          background-color: grey;
          color: white;
          padding: 2px 5px;
        }
        table.summary td {
          text-align: right;
          padding: 2px 5px;
          border-bottom: solid 1px lightgray;
        }
        """))

    # Outputs to include in ZIP file
    zip_contents = []

    # Report one summary section per run
    for run in runs:
        if run is None:
            title = "Summary"
        else:
            title = f"Summary: {run}"
        summary = report.add_section(title)

        # Build the summary based on the discovered outputs
        summary_tbl = Table(("flowcell", "barcode"),
                            flowcell="Flow cell ID",
                            barcode="Barcode ID")
        summary_tbl.add_css_classes("summary")
        if "read_counts" in qc_outputs:
            summary_tbl.append_columns("nreads", nreads="#Reads")
        if "nanoplot" in qc_outputs:
            summary_tbl.append_columns("nanoplot", nanoplot="NanoPlot")
        if "fastq_screen" in qc_outputs:
            summary_tbl.append_columns("screens", screens="Screens")
        summary.add(summary_tbl)

        # Populate the summary table
        previous_fc = None
        for name in fc_bc_names:
            try:
                fc, bc = name.split("_")
            except ValueError:
                fc = name
                bc = "-"
            if fc not in runs[run]:
                continue
            summary_data = {}
            if fc == previous_fc:
                fc = "&nbsp;"
            else:
                previous_fc = fc
            summary_data["flowcell"] = fc
            summary_data["barcode"] = bc
            if "read_counts" in qc_outputs:
                # Add read count
                try:
                    nreads = qc_outputs["read_counts"][name]
                    summary_data["nreads"] = str(nreads)
                except KeyError:
                    pass
            if "nanoplot" in qc_outputs:
                # Add NanoPlot report
                try:
                    nanoplot_report = os.path.join(
                        qc_outputs["nanoplot"][name],
                        "NanoPlot-report.html")
                    summary_data["nanoplot"] = Link(
                        f"NanoPlot report: {name}",
                        os.path.relpath(nanoplot_report, out_dir))
                    zip_contents.append(nanoplot_report)
                except KeyError:
                    pass
            if "fastq_screen" in qc_outputs:
                # Add Fastq Screen output
                try:
                    screen_txt = qc_outputs["fastq_screen"][name]
                    screen_png = os.path.join(
                        "fastq_screen",
                        ".".join(screen_txt.split(".")[:-1]) + ".png")
                    uplot = uscreenplot(
                        [os.path.join(screens_dir, screen_txt)],
                        inline=True)
                    summary_data["screens"] = Img(uplot, href=screen_png)
                    zip_contents.append(os.path.join(out_dir, screen_png))
                except KeyError:
                    pass

            # Add to table
            summary_tbl.add_row(**summary_data)

    # Output to HTML
    html_report = os.path.join(out_dir, "nanoporeqc_report.html")
    report.write(html_report)
    zip_contents.append(html_report)

    # Create a ZIP archive of the report
    if project_name:
        zip_name = f"nanoporeqc-{project_name}.zip"
    else:
        zip_name = f"nanoporeqc.zip"
    zip_path = os.path.join(out_dir, zip_name)
    base_dir = os.path.dirname(out_dir)
    ZipMaker(base_dir, zip_contents).make_archive(zip_path)