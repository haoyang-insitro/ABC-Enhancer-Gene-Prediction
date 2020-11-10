from os.path import basename, join

import numpy as np
import pandas as pd
from redun import Dir, File, script, task
from tools import *

from insitro_core.utils.storage import *


@task()
def make_candidate_regions_from_summits(
    count_file: str,
    macs_peaks: str,
    chrom_sizes: str,
    output_dir: str,
    tmpdir: str,
    regions_whitelist: str = None,
    regions_blacklist: str = None,
    n_enhancers: int = 175000,
    peak_extend: int = 250,
):
    ## Generate enhancer regions from MACS summits
    # 1. Count reads in dhs peaks
    # 2. Take top N regions, get summits, extend summits, merge
    makedirs(output_dir)
    makedirs(tmpdir)
    outfile = join(output_dir, basename(macs_peaks) + ".candidateRegions.bed")
    chrom_sizes_bed = ".".join([chrom_sizes, "bed"])
    # get tmpdir local files
    local_macs_peaks = join(tmpdir, basename(macs_peaks))
    local_chrom_sizes = join(tmpdir, basename(chrom_sizes))
    local_chrom_sizes_bed = ".".join([local_chrom_sizes, "bed"])
    local_regions_whitelist = join(tmpdir, basename(regions_whitelist))
    local_regions_blacklist = join(tmpdir, basename(regions_blacklist))
    local_outfile = join(tmpdir, basename(outfile))
    local_count_file = join(tmpdir, basename(count_file))
    if regions_whitelist:
        whitelist_command = (
            "( bedtools intersect -a "
            + local_regions_whitelist
            + " -b "
            + local_chrom_sizes_bed
            + " -wa | cut -f1-3 && cat ) | "
        )
    else:
        whitelist_command = ""

    if regions_blacklist:
        blacklist_command = (
            "bedtools intersect -v -wa -a stdin -b " + local_regions_blacklist + " | "
        )
    else:
        blacklist_command = ""

    # 2. Take top N regions, get summits, extend summits, merge, remove blacklist, add whitelist, sort and merge
    # use -sorted in intersect command? Not worth it, both files are small
    return script(
        f"""
    #!/bin/bash
    bedtools sort -i {local_count_file} -faidx {local_chrom_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} | \
    bedtools intersect -b stdin -a {local_macs_peaks} -wa | \
    awk '{{print $1"\t"$2 + $10"\t"$2 + $10}}' | \
    bedtools slop -i stdin -b {peak_extend} -g {local_chrom_sizes} | \
    bedtools sort -i stdin -faidx {local_chrom_sizes} | \
    bedtools merge -i stdin |  \
    {blacklist_command} \
    cut -f 1-3 | {whitelist_command}  \
    bedtools sort -i stdin -faidx {local_chrom_sizes} | bedtools merge -i stdin > {local_outfile}
    """,
        inputs=[
            File(count_file).stage(File(local_count_file)),
            File(chrom_sizes).stage(File(local_chrom_sizes)),
            File(chrom_sizes_bed).stage(File(local_chrom_sizes_bed)),
            File(macs_peaks).stage(File(local_macs_peaks)),
            File(regions_whitelist).stage(File(local_regions_whitelist)),
            File(regions_blacklist).stage(File(local_regions_blacklist)),
        ],
        outputs={
            "candidate_enhancer_regions_file": File(outfile).stage(local_outfile),
            "candidate_enhancer_regions_path": outfile,
        },
    )


@task()
def make_candidate_regions_from_peaks(
    count_file: str,
    macs_peaks: str,
    chrom_sizes: str,
    output_dir: str,
    tmpdir: str,
    n_enhancers: int = 175000,
    regions_whitelist: str = None,
    regions_blacklist: str = None,
    peak_extend: int = 250,
    minPeakWidth: int = 500,
):
    makedirs(output_dir)
    makedirs(tmpdir)
    outfile = join(output_dir, basename(macs_peaks) + ".candidateRegions.bed")

    # get tmpdir local files
    local_macs_peaks = join(tmpdir, basename(macs_peaks))
    local_chrom_sizes = join(tmpdir, basename(chrom_sizes))
    local_regions_whitelist = join(tmpdir, basename(regions_whitelist))
    local_regions_blacklist = join(tmpdir, basename(regions_blacklist))
    local_outfile = join(tmpdir, basename(outfile))
    local_count_file = join(tmpdir, basename(count_file))

    ## Generate enhancer regions from MACS narrowPeak - do not use summits
    if regions_whitelist:
        whitelist_command = (
            "(bedtools instersect -a "
            + local_regions_whitelist
            + " -b "
            + local_chrom_sizes_bed
            + " -wa | cut -f1-3 && cat ) | "
        )
    else:
        whitelist_command = ""

    if regions_blacklist:
        blacklist_command = (
            "bedtools intersect -v -wa -a stdin -b " + local_regions_blacklist + " | "
        )
    else:
        blacklist_command = ""

    # 2. Take top N regions, extend peaks (min size 500), merge, remove blacklist, add whitelist, sort and merge
    # use -sorted in intersect command? Not worth it, both files are small
    return script(
        f"""
    #!/bin/bash
    bedtools sort -i {local_count_file} -faidx {local_chrom_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} | \
    bedtools intersect -b stdin -a {local_macs_peaks} -wa | \
    bedtools slop -i stdin -b {peak_extend} -g {local_chrom_sizes} | \
    awk '{{ l=$3-$2; if (l < {minPeakWidth}) {{ $2 = $2 - int(({minPeakWidth}-l)/2); $3 = $3 + int(({minPeakWidth}-l)/2) }} print $1"\t"$2"\t"$3}}' | \
    bedtools sort -i stdin -faidx {local_chrom_sizes} | \
    bedtools merge -i stdin | \
    blacklist_command \
    cut -f 1-3 | {local_whitelist_command} \
    bedtools sort -i stdin -faidx {local_chrom_sizes} | bedtools merge -i stdin > {local_outfile}
    """,
        inputs=[
            File(counts_file).stage(File(local_count_file)),
            File(chrom_sizes).stage(File(local_chrom_sizes)),
            File(macs_peaks).stage(File(local_macs_peaks)),
            File(regions_whitelist).stage(File(local_regions_whitelist)),
            File(regions_blacklist).stage(File(local_regions_blacklist)),
        ],
        outputs={
            "candidate_enhancer_regions_file": File(outfile).stage(local_outfile),
            "candidate_enhancer_regions_path": outfile,
        },
    )
