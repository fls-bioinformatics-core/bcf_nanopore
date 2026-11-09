#!/usr/bin/env python3

# Utility functions for Nanopore QC


import os


def split_fastq_name(fq):
    """
    Split up a Nanopore Fastq filename

    Arguments:
      fq (str): path or filename

    Returns:
      Tuple: tuple of (flowcell, barcode, number)
    """
    fq_name = os.path.basename(fq)
    flowcell_id = fq_name.split("_")[0]
    if "_barcode" in fq_name:
        barcode_id = "barcode" + fq_name.split("_barcode")[1].split("_")[0]
    elif "_unclassified" in fq_name:
        barcode_id = "unclassified"
    else:
        barcode_id = None
    number = fq_name.split(".")[0].split("_")[-1]
    if not number.isdigit():
        number = None
    return (flowcell_id, barcode_id, number)


def group_fastqs(fastqs):
    """
    Group Fastqs by name by flowcell and barcode

    Arguments:
      fastqs (list): list of Fastq names
    """
    flowcells = {}
    for fq in sorted(fastqs):
        flowcell_id, barcode_id, number = split_fastq_name(fq)
        if flowcell_id not in flowcells:
            flowcells[flowcell_id] = {}
        if barcode_id not in flowcells[flowcell_id]:
            flowcells[flowcell_id][barcode_id] = []
        flowcells[flowcell_id][barcode_id].append(fq)
    return flowcells