import os
import traceback
from os.path import basename, join

import numpy as np
import pandas as pd
from neighborhoods import *
from peaks import *
from redun import Dir, File, task
from tools import write_params

from insitro_core.utils.storage import *


@task()
def makeCandidateRegions(
    narrowPeak: str,
    input_bam: str,
    output_dir: str,
    chrom_sizes: str,
    tmpdir: str,
    regions_blacklist: str = None,
    regions_whitelist: str = None,
    peakExtendFromSummit: int = 250,
    nStrongestPeaks: int = 175000,
    ignoreSummits: bool = False,
    minPeakWidth: int = 500,
):
    """
    Inputs:
    narrowPeak: narrowPeak file output by macs2. Must include summits (--call-summits)
    input_bam: DNAase-Seq or atac-Seq input_bam file
    chrom_sizes: File listing chromosome size annotations
    output_dir: output folder where results will be stored; this is created if it doesn't exist
    nStrongestPeaks: Number of peaks to use for defining candidate regions
    peakExtendFromSummit: Number of base pairs to extend each preak from its summit (or from both ends of region if using --ignoreSummits)
    ignoreSummits: Compute peaks using the full peak regions, rather than extending from summit.
    minPeakWidth: Candidate regions whose width is below this threshold are expanded to this width. Only used with --ignoreSummits
    regions_whitelist: Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blacklist
    regions_blacklist: Bed file of regions to forcibly exclude from candidate enhancers
    """
    # create output directory.
    # if output directory is in s3, create local directory from it's basename to serve as workdir
    makedirs(output_dir)
    makedirs(tmpdir)

    write_params(
        {
            "narrowPeak": narrowPeak,
            "input_bam": input_bam,
            "output_dir": output_dir,
            "tmpdir": tmpdir,
            "chrom_sizes": chrom_sizes,
            "regions_blacklist": regions_blacklist,
            "regions_whitelist": regions_whitelist,
            "peakExtendFromSummit": peakExtendFromSummit,
            "nStrongestPeaks": nStrongestPeaks,
            "ignoreSummits": ignoreSummits,
            "minPeakWidth": minPeakWidth,
        },
        output_dir,
        tmpdir,
        "parameters.txt",
    )
    # 1. Count dhs/atac reads in candidate regions
    raw_counts_outfile = join(
        output_dir, basename(narrowPeak) + "." + basename(input_bam) + ".Counts.bed"
    )

    run_count_reads_out = run_count_reads(
        target=input_bam,
        output=raw_counts_outfile,
        output_dir=output_dir,
        tmpdir=tmpdir,
        bed_file=narrowPeak,
        chrom_sizes=chrom_sizes,
        use_fast_count=True,
    )

    # Make candidate regions
    if not ignoreSummits:
        return make_candidate_regions_from_summits(
            count_file=run_count_reads_out["path"],
            macs_peaks=narrowPeak,
            chrom_sizes=chrom_sizes,
            regions_whitelist=regions_whitelist,
            regions_blacklist=regions_blacklist,
            n_enhancers=nStrongestPeaks,
            peak_extend=peakExtendFromSummit,
            output_dir=output_dir,
            tmpdir=tmpdir,
        )
    else:
        return make_candidate_regions_from_peaks(
            count_file=run_count_reads_out["path"],
            macs_peaks=narrowPeak,
            chrom_sizes=chrom_sizes,
            regions_whitelist=regions_whitelist,
            regions_blacklist=regions_blacklist,
            n_enhancers=nStrongestPeaks,
            peak_extend=peakExtendFromSummit,
            minPeakWidth=minPeakWidth,
            output_dir=output_dir,
            tmpdir=tmpdir,
        )
