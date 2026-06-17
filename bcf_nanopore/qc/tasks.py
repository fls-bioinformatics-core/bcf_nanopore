#!/usr/bin/env python3

# Tasks for Nanopore QC pipeline


import os
from .reporting import report_nanopore_qc
from .utils import group_fastqs
from .utils import split_fastq_name
from auto_process_ngs.pipeliner import PipelineTask
from auto_process_ngs.pipeliner import PipelineFunctionTask
from auto_process_ngs.pipeliner import PipelineError
from auto_process_ngs.pipeliner import ListParam
from auto_process_ngs.pipeliner import PipelineParam as Param


class MakeDirs(PipelineTask):
    """
    Set up the output directories

    Arguments:
      out_dir (str): top-level output directory
    """
    def init(self, out_dir):
        self.add_output("out_dir", Param(type='str'))
        self.add_output("merged_fastqs", Param(type='str'))
        self.add_output("nanoplot", Param(type='str'))
        self.add_output("fastq_screen", Param(type='str'))
        self.add_output("read_counts", Param(type='str'))

    def setup(self):
        # Create output directories
        out_dir = self.args.out_dir
        sub_dirs = {d: os.path.join(out_dir, d)
                    for d in ["merged_fastqs",
                              "nanoplot",
                              "fastq_screen",
                              "read_counts"]}
        for name in sub_dirs:
            d = sub_dirs[name]
            if not os.path.exists(d):
                print(f"...making directory '{d}'")
                os.makedirs(d)
        # Set the outputs
        self.output.out_dir.set(out_dir)
        self.output.merged_fastqs.set(sub_dirs["merged_fastqs"])
        self.output.nanoplot.set(sub_dirs["nanoplot"])
        self.output.fastq_screen.set(sub_dirs["fastq_screen"])
        self.output.read_counts.set(sub_dirs["read_counts"])


class GetFastqs(PipelineTask):
    """
    Get the Fastqs to run the QC on
    """

    def init(self, fastqs_in):
        self.add_output("fastqs", ListParam())

    def setup(self):
        # Initially just copy the input list
        for fq in self.args.fastqs_in:
            self.output.fastqs.append(fq)


class MergeFastqs(PipelineTask):
    """
    Merge subset of Fastqs within flow cell and barcode groups
    """

    def init(self, fastqs, out_dir):
        self.conda("seqtk=1.5")
        self.add_output("fastqs", ListParam())

    def setup(self):
        fq_groups = group_fastqs(self.args.fastqs)
        if self.args.fastqs[0].endswith(".fastq.gz"):
            ftype = "fastqgz"
        elif self.args.fastqs[0].endswith(".fastq"):
            ftype = "fastq"
        else:
            raise Exception("Unknown extension for Fastqs")
        for flowcell in fq_groups:
            barcodes = fq_groups[flowcell]
            for barcode in barcodes:
                fqs = barcodes[barcode]
                if barcode is None:
                    name = f"{flowcell}"
                else:
                    name = f"{flowcell}_{barcode}"
                print(f"-- {name}: {len(fqs)} Fastqs")
                fq_basename = "_".join(os.path.basename(fqs[0]).
                                       split(".")[0].
                                       split("_")[:-1]) + ".fastq.gz"
                fastq_out = os.path.join(self.args.out_dir, fq_basename)
                if not os.path.exists(fastq_out):
                    print(f"...creating '{fq_basename}'")
                    self.add_cmd(
                        f"Merge random Fastq subset for '{name}'",
                        """
                        # Make temporary working directory
                        WORKDIR=$(mktemp -d --tmpdir=. {name}.XXX)
                        cd $WORKDIR
                        # Get subset of reads from each Fastq
                        fastqs="{fastqs}"
                        idx=0
                        for fq in $fastqs ; do
                           idx=$((idx+1))
                           fq_out=subset_${{idx}}.fq
                           seqtk sample -s 1234 $fq {subset} >$fq_out
                        done
                        # Merge Fastqs
                        cat subset_*.fq | pigz >merged.fastq.gz
                        # Copy to final location
                        /bin/cp merged.fastq.gz {fastq_out}
                        """.format(name=name,
                                   subset=1000,
                                   fastqs=" ".join(fqs),
                                   fastq_out=fastq_out))
                else:
                    print(f"... found '{fq_basename}'")

    def finish(self):
        for fq in os.listdir(self.args.out_dir):
            self.output.fastqs.append(
                os.path.join(self.args.out_dir, fq))


class GetReadCounts(PipelineTask):
    """
    Get read counts for groups of Fastqs
    """

    def init(self, fastqs, out_dir):
        # Location of outputs
        self.add_output("out_dir", Param(type="str"))
        self.add_output("counts", Param())

    def setup(self):
        fq_groups = group_fastqs(self.args.fastqs)
        if self.args.fastqs[0].endswith(".fastq.gz"):
            ftype = "fastqgz"
        elif self.args.fastqs[0].endswith(".fastq"):
            ftype = "fastq"
        else:
            raise Exception("Unknown extension for Fastqs")
        for fc in fq_groups:
            barcodes = fq_groups[fc]
            for bc in barcodes:
                if bc is not None:
                    name = f"{fc}_{bc}"
                else:
                    name = fc
                count_file = os.path.join(self.args.out_dir,
                                          f"{name}_count.out")
                if os.path.exists(count_file):
                    # Already have the counts
                    continue
                self.add_cmd(
                    # Use unpigz instead of zcat for compressed Fastqs
                    # See https://unix.stackexchange.com/a/363739
                    f"Read counts for {name}",
                    """
                    echo $(($({cmd} {fastqs} | wc -l)/4)) >{count_file}
                    """.format(cmd=("unpigz -c" if ftype == "fastqgz"
                                                        else "cat"),
                               fastqs=" ".join(barcodes[bc]),
                               count_file=count_file))

    def finish(self):
        fq_groups = group_fastqs(self.args.fastqs)
        counts = {}
        for fc in fq_groups:
            barcodes = fq_groups[fc]
            for bc in barcodes:
                if bc is not None:
                    name = f"{fc}_{bc}"
                else:
                    name = f"{fc}"
                count_file = os.path.join(self.args.out_dir,
                                          f"{name}_count.out")
                if os.path.exists(count_file):
                    with open(count_file, "rt") as log:
                        nreads = int(log.read().strip())
                        print(f"{name}: {nreads}")
                        counts[name] = nreads
        self.output.out_dir.set(self.args.out_dir)
        self.output.counts.set(counts)


class FilterFastqs(PipelineTask):
    """
    Filter Fastq list by read count
    """

    def init(self, fastqs, counts, min_reads=1000):
        self.add_output("fastqs", ListParam())

    def setup(self):
        for fq in self.args.fastqs:
            fc, bc, number = split_fastq_name(fq)
            if bc is not None:
                name = f"{fc}_{bc}"
            else:
                name = fc
            if self.args.counts[name] >= self.args.min_reads:
                self.output.fastqs.append(fq)


class FastqScreen(PipelineTask):
    """
    Run FastqScreen

    Cannabilised from auto-process-ngs/qc/modules/fastq_screen

    Arguments:
      fastqs (list): list of paths to Fastq files to
        run FastqScreen on
      out_dir (str): directory for outputs
      conf_file (str): path to conf file for FastqScreen
      screen_name (str): name for screen (defaults to the
        name of the conf file)
      subset (int): explicitly specify the subset size
        for running Fastq_screen
      nthreads (int): number of threads/processors to
        use (defaults to number of slots set in runner)
      aligner (str): aligner to use, either "bowtie" or
        "minimap2" (default is "minimap2")
    """

    def init(self, fastqs, out_dir, conf_file, screen_name=None,
             subset=None, nthreads=None, aligner="minimap2"):
        self.conda("fastq-screen=0.16.0",
                   "perl-gd")
        if aligner == "bowtie":
            # Also need to specify tbb=2020.2 for bowtie
            # See https://www.biostars.org/p/494922/
            self.conda("bowtie=1.2.3",
                       "tbb=2020.2")
        elif aligner == "minimap2":
            self.conda("minimap2=2.30")
        else:
            raise PipelineError(f"'{aligner}': unrecognised aligner")
        # Location of outputs
        self.add_output("out_dir", Param(type="str"))

    def setup(self):
        # Fastqs
        print(f"Fastqs: {self.args.fastqs}")
        if not self.args.fastqs:
            print("Nothing to do")
            return
        else:
            fastqs = sorted(self.args.fastqs)
        # Conf file
        print(f"Conf file = {self.args.conf_file}")
        if not os.path.isfile(self.args.conf_file):
            raise Exception(f"{self.args.conf_file}: conf file not "
                            "found (or is not a file)")
        # Screen name
        if self.args.screen_name:
            screen_name = str(self.args.screen_name)
        else:
            screen_name = os.path.basename(self.args.conf_file).split(".")[0]
        # Set up the FastqScreen runs for each Fastq
        for fastq in fastqs:
            # Base name for Fastq file
            fastq_basename = os.path.basename(fastq)
            while fastq_basename.split('.')[-1] in ('fastq', 'fq', 'gz'):
                fastq_basename = '.'.join(fastq_basename.split('.')[:-1])
            # Determine base names for outputs
            screen_basename = "%s_screen_%s" % (fastq_basename,
                                                screen_name)
            # Check if outputs already exist
            outputs_exist = True
            for ext in ("png", "txt"):
                if not os.path.exists(os.path.join(
                        self.args.out_dir,
                        f"{screen_basename}.{ext}")):
                    outputs_exist = False
                    break
            if outputs_exist:
                print(f"...skipping {os.path.basename(fastq)}")
                continue
            # Set parameters
            if self.args.nthreads:
                nthreads = self.args.nthreads
            else:
                nthreads = self.runner_nslots
            if self.args.subset is not None:
                subset_option = "--subset %d" % self.args.subset
            else:
                subset_option = ""
            # Run FastqScreen
            self.add_cmd(
                "Run FastqScreen '%s' on %s" % (screen_name,
                                                os.path.basename(fastq)),
                """
                # Make temporary working directory
                WORKDIR=$(mktemp -d --tmpdir=. {name}.XXX)
                # Run FastqScreen: {screen_name}
                fastq_screen \\
                    --aligner {aligner} \\
                    --conf {fastq_screen_conf} \\
                    --threads {nthreads} {subset_option} \\
                    --outdir $WORKDIR --force \\
                    {fastq}
                # Rename and move outputs to final location
                out_txt=$WORKDIR/{fastq_basename}_screen.txt
                if [ -e $out_txt ] ; then
                    /bin/mv $out_txt {out_dir}/{screen_basename}.txt
                else
                    echo "ERROR missing $out_txt" >&2
                    exit 1
                fi
                out_png=$WORKDIR/{fastq_basename}_screen.png
                if [ -e $out_png ] ; then
                    /bin/mv $out_png {out_dir}/{screen_basename}.png
                else
                    echo "ERROR missing $out_png output" >&2
                    exit 1
                fi
                """.format(name=screen_basename,
                           screen_name=screen_name,
                           fastq=fastq,
                           aligner=self.args.aligner,
                           out_dir=self.args.out_dir,
                           fastq_screen_conf=self.args.conf_file,
                           nthreads=nthreads,
                           subset_option=subset_option,
                           fastq_basename=fastq_basename,
                           screen_basename=screen_basename))

    def finish(self):
        # Set output
        self.output.out_dir.set(self.args.out_dir)


class NanoPlot(PipelineTask):
    """
    Run NanoPlot

    Arguments:
      fastqs (list): input Fastqs
      out_dir (str): path to output directory
      nthreads (int): number of threads to run NanoPlot with
        (default: 1)
    """

    def init(self, fastqs, out_dir, nthreads=1):
        self.conda("nanoplot=1.46.2")
        self.add_output("out_dir", Param(type="str"))

    def setup(self):
        # Expected outputs
        nanoplot_outputs = [
            "NanoPlot-report.html",
            "NanoStats.txt",
            "NanoPlot-data.tsv.gz",
            "*.png"
        ]
        # Group Fastqs
        fq_groups = group_fastqs(self.args.fastqs)
        # Loop over flowcell/barcode combinations
        for flowcell in fq_groups:
            barcodes = fq_groups[flowcell]
            for barcode in barcodes:
                # Fastqs in this group
                fqs = sorted(barcodes[barcode])
                # Create group name
                if barcode is None:
                    group_name = f"{flowcell}"
                else:
                    group_name = f"{flowcell}_{barcode}"
                # Check if Nanoplot outputs already exist for the group
                out_dir = os.path.join(self.args.out_dir, group_name)
                run_nanoplot = False
                for output in [os.path.join(out_dir, f)
                               for f in nanoplot_outputs
                               if "*" not in f]:
                    if not os.path.exists(output):
                        print(f"...missing Nanoplot outputs for "
                              f"{group_name}")
                        run_nanoplot = True
                        break
                if not run_nanoplot:
                    print(f"...outputs already present for {group_name}")
                    continue
                # Run NanoPlot
                self.add_cmd(
                    f"Run NanoPlot for '{group_name}'",
                    """
                    # Make temporary working directory
                    WORKDIR=$(mktemp -d --tmpdir=. {name}.XXX)
                    cd $WORKDIR
                    # Run NanoPlot
                    NanoPlot \\
                        -t {threads} \\
                        --fastq {fastqs} \\
                        --raw
                    # Make output directory
                    if [ ! -e {out_dir} ] ; then
                        mkdir -p {out_dir}
                    fi
                    # Copy outputs
                    OUTPUTS="{nanoplot_outputs}"
                    for output in $OUTPUTS ; do
                        /bin/cp $output {out_dir}
                    done
                    """.format(name=group_name,
                               threads=self.args.nthreads,
                               fastqs=" ".join(fqs),
                               nanoplot_outputs=" ".join(nanoplot_outputs),
                               out_dir=out_dir))

    def finish(self):
        # Set output
        self.output.out_dir.set(self.args.out_dir)


class Report(PipelineFunctionTask):
    """
    Create a HTML report

    (NB currently the fastq_screen, nanoplot etc arguments
    are "flag" parameters used to check if the required
    tasks have completed.)

    Arguments:
      out_dir (str): top-level output directory
      project_dir (str): PromethION project directory (optional)
      fastq_screen (str): path to FastqScreen outputs
      nanoplot (str): path to NanoPlot outputs
      read_counts (str): path to read count outputs
    """

    def init(self, out_dir, project_dir=None, fastq_screen=None,
             nanoplot=None, read_counts=None):
        pass

    def setup(self):
        # Generate the QC report
        self.add_call("Generate HTML report",
                      report_nanopore_qc,
                      self.args.out_dir,
                      project_dir=self.args.project_dir)