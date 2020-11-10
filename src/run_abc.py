import os
import sys
from os.path import basename

from redun import Dir, File, Scheduler, namespace, task

# import biotools & utils scripts one directory up
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "biotools_redun"))
from biotools import *
from makeCandidateRegions import *
from neighborhoods import *
from predict import *


@task()
def workflow(
    output_dir: str,
    candidate_enhancer_regions: str,
    chrom_sizes: str,
    enhancers: str,
    genes_file: str,
    hic_resolution: int,
    input_bam: str,
    tmpdir: str,
    atac: str = "",
    bucket_region: str = "us-west-2",
    celltype: str = None,
    chromosomes: str = "all",
    default_accessibility_feature: str = None,
    dhs: str = "",
    enhancer_class_override: str = None,
    expression_cutoff: float = 1,
    expression_table: str = "",
    gene_name_annotations: str = "symbol",
    genes_for_class_assignment: str = None,
    genome: str = "hs",
    h3k27ac: str = "",
    hicdir: str = "",
    hic_gamma: float = 0.87,
    hic_gamma_reference: float = 0.87,
    hic_pseudocount_distance: int = 1e6,
    hic_type: str = "juicebox",
    hic_is_doubly_stochastic: bool = False,
    ignoreSummits: bool = False,
    include_chrY: bool = True,
    make_all_putative: bool = False,
    minPeakWidth: int = 500,
    nStrongestPeaks: int = 175000,
    narrowPeak: str = None,
    peakExtendFromSummit: int = 250,
    primary_gene_identifier: str = "symbol",
    promoter_activity_quantile_cutoff: float = 0.4,
    pthresh: float = 0.1,
    qnorm: str = None,
    regions_blacklist: str = None,
    regions_whitelist: str = None,
    run_all_genes: bool = False,
    scale_hic_using_powerlaw: bool = False,
    score_column: str = "ABC.Score",
    skip_gene_counts: bool = False,
    skip_rpkm_quantile: bool = False,
    supplementary_features: str = None,
    threshold: float = 0.022,
    tss_hic_contribution: float = 100,
    tss_slop: int = 500,
    ubiquitously_expressed_genes: str = None,
    use_secondary_counting_method: bool = True,
    use_hdf5: bool = False,
    window=5e6,
):
    """
    Required Arguments
    ==================

    makeCandidateRegions
    ===================
    narrowPeak:str=None
        narrowPeak file output by macs2. Must include summits (--call-summits)
    genes_file:str
        bed file with gene annotations. Must be in bed-6 format. Will be used to assign TSS to genes.
    ubiquitously_expressed_genes:str

    """
    if narrowPeak is None:
        # perform peak calling on input bam with MACS2
        print("running macs2")
        run_name = basename(input_bam).rstrip(".bam")
        macs2_peaks_from_bam_out = macs2_peaks_from_bam(
            input_bam=input_bam,
            pthresh=pthresh,
            genome=genome,
            output_dir=output_dir,
            tmpdir=tmpdir,
            run_name=run_name,
            bucket_region=bucket_region,
        )

        # sort the bed file
        print("sorting the bed file")
        bedtools_sort_out = bedtools_sort(
            macs2_peaks_from_bam_out["narrowPeak"],
            tmpdir=tmpdir,
            bucket_region=bucket_region,
        )
    else:
        bedtools_sort_out = bedtools_sort(
            narrowPeak,
            tmpdir=tmpdir,
            bucket_region=bucket_region,
        )
    # run makeCandidateRegions task
    print("making candidate regions")
    makeCandidateRegions_out = makeCandidateRegions(
        narrowPeak=bedtools_sort_out["path"],
        input_bam=input_bam,
        output_dir=output_dir,
        tmpdir=tmpdir,
        chrom_sizes=chrom_sizes,
        regions_blacklist=regions_blacklist,
        regions_whitelist=regions_whitelist,
        peakExtendFromSummit=peakExtendFromSummit,
        nStrongestPeaks=nStrongestPeaks,
        ignoreSummits=ignoreSummits,
        minPeakWidth=minPeakWidth,
    )
    print("running processCellType")
    processCellType_out = processCellType(
        genes_file=genes_file,
        ue_file=ubiquitously_expressed_genes,
        chrom_sizes=chrom_sizes,
        output_dir=output_dir,
        tmpdir=tmpdir,
        default_accessibility_feature=default_accessibility_feature,
        expression_table=expression_table,
        gene_name_annotations=gene_name_annotations,
        primary_gene_identifier=primary_gene_identifier,
        celltype=celltype,
        atac=atac,
        dhs=dhs,
        h3k27ac=h3k27ac,
        supplementary_features=supplementary_features,
        class_gene_file=genes_for_class_assignment,
        candidate_enhancer_regions=makeCandidateRegions_out[
            "candidate_enhancer_regions_path"
        ],
        skip_gene_counts=skip_gene_counts,
        use_secondary_counting_method=use_secondary_counting_method,
        skip_rpkm_quantile=skip_rpkm_quantile,
        qnorm=qnorm,
        tss_slop=tss_slop,
        enhancer_class_override=enhancer_class_override,
    )

    print("running prediction")
    return predict(
        output_dir=output_dir,
        tmpdir=tmpdir,
        enhancers=processCellType_out["enhancers"][0],
        genes=processCellType_out["genes_for_class_assignment"],
        celltype=celltype,
        hic_resolution=hic_resolution,
        window=window,
        score_column=score_column,
        threshold=threshold,
        hicdir=hicdir,
        tss_hic_contribution=tss_hic_contribution,
        hic_pseudocount_distance=hic_pseudocount_distance,
        hic_type=hic_type,
        hic_is_doubly_stochastic=hic_is_doubly_stochastic,
        scale_hic_using_powerlaw=scale_hic_using_powerlaw,
        hic_gamma=hic_gamma,
        hic_gamma_reference=hic_gamma_reference,
        run_all_genes=run_all_genes,
        expression_cutoff=expression_cutoff,
        promoter_activity_quantile_cutoff=promoter_activity_quantile_cutoff,
        make_all_putative=make_all_putative,
        use_hdf5=use_hdf5,
        tss_slop=tss_slop,
        chromosomes=chromosomes,
        include_chrY=include_chrY,
    )
