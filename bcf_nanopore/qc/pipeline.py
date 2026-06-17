#!/usr/bin/env python3

# Create pipeline for running FastqScreen on Nanopore data


import os
import shutil
import tempfile
from .tasks import MakeDirs
from .tasks import GetFastqs
from .tasks import GetReadCounts
from .tasks import FilterFastqs
from .tasks import MergeFastqs
from .tasks import FastqScreen
from .tasks import NanoPlot
from .tasks import Report
from bcftbx.JobRunner import SimpleJobRunner
from auto_process_ngs.pipeliner import Pipeline
from auto_process_ngs.pipeliner import PipelineFailure


class NanoporeQC(Pipeline):
    """
    Pipeline for running QC on Nanopore data

    Arguments:
        modules (list): list of QC modules to execute
    """
    def __init__(self, modules=["nanoplot", "fastq_screen"]):
        # Initialise the pipeline superclass
        Pipeline.__init__(self, name="NanoporeQC")
        # Define parameters
        self.add_param("fastqs", type=list)
        self.add_param("nthreads", type=int)
        self.add_param("fastq_subset", type=int)
        self.add_param("fastq_screen_conf_file", type=str)
        self.add_param("read_counts_out_dir", type=str)
        self.add_param("fastq_screen_out_dir", type=str)
        self.add_param("nanoplot_out_dir", type=str)
        self.add_param("project_dir", type=str)
        self.add_param("out_dir", type=str)

        # Define runners
        self.add_runner("fastq_screen")
        self.add_runner("nanoplot")

        # Add tasks to the pipeline
        make_dirs = MakeDirs(
            "Create output directories",
            self.params.out_dir)
        self.add_task(make_dirs)

        get_fastqs = GetFastqs(
            "Fetch Fastqs",
            self.params.fastqs)
        self.add_task(get_fastqs)

        get_read_counts = GetReadCounts(
            "Get read counts by flowcell and barcode",
            get_fastqs.output.fastqs,
            make_dirs.output.read_counts)
        self.add_task(get_read_counts)
        self.params["read_counts_out_dir"] = \
            get_read_counts.output.out_dir

        filter_fastqs = FilterFastqs(
            "Filter Fastqs on number of reads",
            get_fastqs.output.fastqs,
            get_read_counts.output.counts)
        self.add_task(filter_fastqs)

        merge_fastqs = MergeFastqs(
            "Merge Fastqs across flow cells and barcodes",
            filter_fastqs.output.fastqs,
            make_dirs.output.merged_fastqs)
        self.add_task(merge_fastqs)

        if "fastq_screen" in modules:
            run_fastq_screen = FastqScreen(
                "Run FastqScreen",
                merge_fastqs.output.fastqs,
                make_dirs.output.fastq_screen,
                self.params.fastq_screen_conf_file,
                subset=self.params.fastq_subset,
                nthreads=self.params.nthreads,
                aligner="minimap2")
            self.add_task(run_fastq_screen,
                          runner=self.runners["fastq_screen"])
            self.params["fastq_screen_out_dir"] = \
                run_fastq_screen.output.out_dir

        if "nanoplot" in modules:
            run_nanoplot = NanoPlot(
                "Run Nanoplot",
                filter_fastqs.output.fastqs,
                make_dirs.output.nanoplot,
                nthreads=self.params.nthreads
            )
            self.add_task(run_nanoplot,
                          runner=self.runners["nanoplot"])
            self.params["nanoplot_out_dir"] = \
                run_nanoplot.output.out_dir

        make_report = Report(
            "Make HTML report",
            self.params.out_dir,
            project_dir=self.params.project_dir,
            fastq_screen=self.params.fastq_screen_out_dir,
            nanoplot=self.params.nanoplot_out_dir,
            read_counts=self.params.read_counts_out_dir)
        self.add_task(make_report)

    def run(self, fastqs, fastq_screen_conf_file, out_dir,
            nthreads=8, project_dir=None):
        """
        Run the pipeline
        """
        working_dir = tempfile.mkdtemp(prefix="__nanoporeqc.",
                                       suffix=".tmp",
                                       dir=os.getcwd())
        status = Pipeline.run(
            self,
            working_dir=working_dir,
            params={
                "project_dir": project_dir,
                "fastqs": fastqs,
                "fastq_screen_conf_file": fastq_screen_conf_file,
                "out_dir": out_dir,
                "nthreads": nthreads,
            },
            max_jobs=12,
            max_slots=12,
            enable_conda=True,
            conda_env_dir=os.path.join(os.getcwd(),
                                       "__nanoporeqc.conda"),
            runners={
                "fastq_screen": SimpleJobRunner(nslots=nthreads),
                "nanoplot": SimpleJobRunner(nslots=nthreads)
            },
            default_runner=SimpleJobRunner(),
            verbose=False,
            exit_on_failure=PipelineFailure.DEFERRED)
        if status == 0:
            shutil.rmtree(working_dir)
        return status
