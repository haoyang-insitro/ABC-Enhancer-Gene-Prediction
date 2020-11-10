import linecache
import os
import pdb
import time
import traceback
from os.path import basename, join

import numpy as np
import pandas as pd
import pyranges as pr
from redun import Dir, File, script, task
from scipy import interpolate
from tools import *

from insitro_core.utils.cloud.bucket_utils import download_file
from insitro_core.utils.storage import *

pd.options.display.max_colwidth = (
    10000  # seems to be necessary for pandas to read long file names... strange
)


@task()
def processCellType(
    genes_file: str,
    ue_file: str,
    chrom_sizes: str,
    output_dir: str,
    tmpdir: str,
    expression_table: str,
    gene_name_annotations: str,
    primary_gene_identifier: str,
    celltype: str,
    atac: str,
    dhs: str,
    h3k27ac: str,
    supplementary_features: str,
    class_gene_file: str,
    default_accessibility_feature: str,
    skip_gene_counts: bool,
    use_secondary_counting_method: bool,
    skip_rpkm_quantile: bool,
    qnorm: bool,
    candidate_enhancer_regions,
    enhancer_class_override: str,
    tss_slop: int = 500,
):
    # run run.neighborhoods task
    print("loading genes")
    if expression_table is not None:
        expression_table_list = expression_table.split(",")
    else:
        expression_table_list = ""

    load_genes_out = load_genes(
        genes_file=genes_file,
        ue_file=ue_file,
        chrom_sizes=chrom_sizes,
        output_dir=output_dir,
        tmpdir=tmpdir,
        expression_table_list=expression_table_list,
        gene_name_annotations=gene_name_annotations,
        primary_gene_identifier=primary_gene_identifier,
        celltype=celltype,
        class_gene_file=class_gene_file,
    )
    genes_for_class_assignment = load_genes_out["genes_for_class_assignment"]
    genes = load_genes_out["genes"]
    default_accessibility_feature = determine_accessibility_feature(
        default_accessibility_feature, atac, dhs
    )
    features = get_features(h3k27ac, atac, dhs, supplementary_features)
    if skip_gene_counts is False:
        print("annotating genes with features")
        genes_for_class_assignment = annotate_genes_with_features(
            genes=load_genes_out["genes"],
            chrom_sizes=chrom_sizes,
            use_fast_count=(not use_secondary_counting_method),
            default_accessibility_feature=default_accessibility_feature,
            features=features,
            output_dir=output_dir,
            tmpdir=tmpdir,
        )

    # Setup Candidate Enhancers
    print("loading enhancers")
    enhancers = load_enhancers(
        genes=load_genes_out["genes_for_class_assignment"],
        chrom_sizes=chrom_sizes,
        candidate_enhancer_regions=candidate_enhancer_regions,
        skip_rpkm_quantile=skip_rpkm_quantile,
        qnorm=qnorm,
        tss_slop_for_class_assignment=tss_slop,
        use_fast_count=(not use_secondary_counting_method),
        default_accessibility_feature=default_accessibility_feature,
        features=features,
        celltype=celltype,
        class_override_file=enhancer_class_override,
        output_dir=output_dir,
        tmpdir=tmpdir,
    )
    print("Neighborhoods Complete! \n")
    return {
        "genes": genes,
        "genes_for_class_assignment": genes_for_class_assignment,
        "enhancers": enhancers,
    }


@task()
def load_genes(
    genes_file: str,
    ue_file: str,
    chrom_sizes: str,
    output_dir: str,
    tmpdir: str,
    expression_table_list: str,
    gene_name_annotations: str = "symbol",
    primary_gene_identifier: str = "symbol",
    celltype: str = None,
    class_gene_file: str = None,
):
    """
    genes_file: bed file with gene annotations. Must be in bed-6 format. Will be used to assign TSS to genes.
    """
    makedirs(tmpdir)
    makedirs(output_dir)
    print("made dirs")
    genes = process_gene_bed(
        bed=genes_file,
        tmpdir=tmpdir,
        name_cols=gene_name_annotations,
        main_name=primary_gene_identifier,
        chrom_sizes=chrom_sizes,
    )
    return process_genes(
        genes=genes,
        tmpdir=tmpdir,
        expression_table_list=expression_table_list,
        primary_gene_identifier=primary_gene_identifier,
        ue_file=ue_file,
        celltype=celltype,
        output_dir=output_dir,
        class_gene_file=class_gene_file,
    )


@task()
def process_genes(
    genes,
    tmpdir,
    expression_table_list,
    primary_gene_identifier,
    ue_file,
    celltype,
    output_dir,
    class_gene_file,
):
    genes[["chr", "start", "end", "name", "score", "strand"]].to_csv(
        join(tmpdir, "GeneList.bed"), sep="\t", index=False, header=False
    )
    # upload genes to outputdir
    upload_file(join(tmpdir, "GeneList.bed"), output_dir, overwrite_ok=True)
    if len(expression_table_list) > 0:
        # Add expression information
        names_list = []
        print("Using gene expression from files: {} \n".format(expression_table_list))

        for expression_table in expression_table_list:
            try:
                download_file(expression_table, tmpdir, overwrite_ok=True)
                name = basename(expression_table)
                local_expression_table = join(tmpdir, basename(expression_table))
                expr = pd.read_table(
                    local_expression_table,
                    names=[primary_gene_identifier, name + ".Expression"],
                )
                expr[name + ".Expression"] = expr[name + ".Expression"].astype(float)
                expr = expr.groupby(primary_gene_identifier).max()

                genes = genes.merge(
                    expr, how="left", right_index=True, left_on="symbol"
                )
                names_list.append(name + ".Expression")
            except Exception as e:
                print(e)
                traceback.print_exc()
                print("Failed on {}".format(expression_table))
        genes["Expression"] = genes[names_list].mean(axis=1)
        genes["Expression.quantile"] = genes["Expression"].rank(
            method="average", na_option="top", ascending=True, pct=True
        )
    else:
        genes["Expression"] = np.NaN

    # Ubiquitously expressed annotation
    if ue_file is not None:
        download_file(ue_file, tmpdir, overwrite_ok=True)
        local_ue_file = join(tmpdir, basename(ue_file))
        ubiq = pd.read_csv(local_ue_file, sep="\t")
        genes["is_ue"] = genes["name"].isin(ubiq.iloc[:, 0].values.tolist())

    # cell type
    genes["celltype"] = celltype

    if class_gene_file is None:
        genes_for_class_assignment = genes
    else:
        genes_for_class_assignment = process_gene_bed(
            genes_for_class_assignment,
            gene_name_annotations,
            primary_gene_identifier,
            chrom_sizes,
            fail_on_nonunique=False,
            tmpdir=tmpdir,
        )

    return {"genes": genes, "genes_for_class_assignment": genes_for_class_assignment}


@task()
def annotate_genes_with_features(
    genes: dict,
    chrom_sizes: str,
    output_dir: str,
    tmpdir: str,
    skip_gene_counts: bool = False,
    features: dict = {},
    force: bool = False,
    use_fast_count: bool = True,
    default_accessibility_feature: str = "",
):
    # pull down a local copy of chromsizes, as that will be used in several locations downstream
    download_file(chrom_sizes, tmpdir, overwrite_ok=True)
    download_file(chrom_sizes + ".bed", tmpdir, overwrite_ok=True)

    chrom_sizes = join(tmpdir, basename(chrom_sizes))
    # Setup files for counting
    bounds_bed = join(output_dir, "GeneList.bed")
    tss1kb = make_tss_region_file(genes, output_dir, tmpdir, chrom_sizes)

    # Count features over genes and promoters
    genes = count_features_for_bed(
        genes,
        bounds_bed,
        chrom_sizes,
        features,
        tmpdir,
        output_dir,
        "Genes",
        force=force,
        use_fast_count=use_fast_count,
    )
    tsscounts = count_features_for_bed(
        tss1kb["df"],
        tss1kb["path"],
        chrom_sizes,
        features,
        tmpdir,
        output_dir,
        "Genes.TSS1kb",
        force=force,
        use_fast_count=use_fast_count,
    )
    return merge_genes_tss(
        genes, tsscounts, tmpdir, output_dir, default_accessibility_feature
    )


@task()
def merge_genes_tss(
    genes, tsscounts, tmpdir, output_dir, default_accessibility_feature
):
    tsscounts = tsscounts.drop(["chr", "start", "end", "score", "strand"], axis=1)
    merged = genes.merge(tsscounts, on="name", suffixes=["", ".TSS1Kb"])

    access_col = default_accessibility_feature + ".RPKM.quantile.TSS1Kb"

    if "h3k27ac.RPKM.quantile.TSS1Kb" in merged.columns:
        merged["PromoterActivityQuantile"] = (
            (0.0001 + merged["h3k27ac.RPKM.quantile.TSS1Kb"])
            * (0.0001 + merged[access_col])
        ).rank(method="average", na_option="top", ascending=True, pct=True)
    else:
        merged["PromoterActivityQuantile"] = ((0.0001 + merged[access_col])).rank(
            method="average", na_option="top", ascending=True, pct=True
        )
    merged.to_csv(
        join(tmpdir, "GeneList.txt"),
        sep="\t",
        index=False,
        header=True,
        float_format="%.6f",
    )
    upload_file(join(tmpdir, "GeneList.txt"), output_dir, overwrite_ok=True)
    return merged


@task()
def make_tss_region_file(genes, output_dir, tmpdir, sizes, tss_slop=500):
    # Given a gene file, define 1kb regions around the tss of each gene
    sizes_pr = df_to_pyranges(
        pd.read_csv(sizes + ".bed", sep="\t", header=None).rename(
            columns={0: "chr", 1: "start", 2: "end"}
        )
    )
    tss1kb = genes.loc[:, ["chr", "start", "end", "name", "score", "strand"]]
    tss1kb["start"] = genes["tss"]
    tss1kb["end"] = genes["tss"]
    tss1kb = df_to_pyranges(tss1kb).slack(tss_slop)
    tss1kb = pr.gf.genome_bounds(tss1kb, sizes_pr).df[
        ["Chromosome", "Start", "End", "name", "score", "strand"]
    ]
    tss1kb.columns = ["chr", "start", "end", "name", "score", "strand"]
    tss1kb.sort_values(["chr", "start", "end"])
    tss1kb_file = os.path.join(tmpdir, "GeneList.TSS1kb.bed")
    tss1kb.to_csv(tss1kb_file, header=False, index=False, sep="\t")

    local_chrom_sizes = join(tmpdir, basename(sizes))
    tss1kb_out_file = join(output_dir, "GeneList.TSS1kb.bed")

    return script(
        f"""
        bedtools sort -faidx {local_chrom_sizes} -i {tss1kb_file} > {tss1kb_file}.sorted; mv {tss1kb_file}.sorted {tss1kb_file}
        """,
        inputs=[File(join(output_dir, basename(sizes))).stage(File(local_chrom_sizes))],
        outputs={
            "file": File(tss1kb_out_file).stage(tss1kb_file),
            "path": tss1kb_out_file,
            "df": tss1kb,
        },
    )


@task()
def process_gene_bed(
    bed,
    tmpdir: str,
    name_cols: str = "symbol",
    main_name: str = "symbol",
    chrom_sizes: str = None,
    fail_on_nonunique: bool = True,
):
    # get local files from s3
    makedirs(tmpdir)
    local_bed = join(tmpdir, basename(bed))
    print(str(local_bed))
    download_file(bed, tmpdir, overwrite_ok=True)

    if chrom_sizes is not None:
        local_chrom_sizes = join(tmpdir, basename(chrom_sizes))
        print(str(local_chrom_sizes))
        download_file(chrom_sizes, tmpdir, overwrite_ok=True)
    else:
        local_chrom_sizes = None
    bed = read_bed(local_bed)
    return process_gene_bed_filter(
        bed, name_cols, main_name, local_chrom_sizes, fail_on_nonunique
    )


@task()
def process_gene_bed_filter(
    bed, name_cols, main_name, local_chrom_sizes, fail_on_nonunique
):
    try:
        bed = bed.drop(
            [
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts",
            ],
            axis=1,
        )
    except Exception as e:
        pass
    assert main_name in name_cols

    names = bed.name.str.split(";", expand=True)
    assert len(names.columns) == len(name_cols.split(","))
    names.columns = name_cols.split(",")
    bed = pd.concat([bed, names], axis=1)

    bed["name"] = bed[main_name]
    # bed = bed.sort_values(by=['chr','start']) #JN Keep original sort order

    bed["tss"] = get_tss_for_bed(bed)

    bed.drop_duplicates(inplace=True)

    # Remove genes that are not defined in chromosomes file
    if local_chrom_sizes is not None:
        sizes = pd.read_csv(local_chrom_sizes, sep="\t", header=None)
        sizes = sizes.rename(columns={0: "chr", 1: "size"})
        bed["chr"] = bed["chr"].astype(
            "str"
        )  # JN needed in case chromosomes are all integer
        bed = bed[bed["chr"].isin(set(sizes["chr"].values))]

    # Enforce that gene names should be unique
    if fail_on_nonunique:
        assert len(set(bed["name"])) == len(
            bed["name"]
        ), "Gene IDs are not unique! Failing. Please ensure unique identifiers are passed to --genes"
    return bed


def get_tss_for_bed(bed):
    assert_bed3(bed)
    tss = bed["start"].copy()
    tss.loc[bed.loc[:, "strand"] == "-"] = bed.loc[bed.loc[:, "strand"] == "-", "end"]

    return tss


def assert_bed3(df):
    assert type(df).__name__ == "DataFrame"
    assert "chr" in df.columns
    assert "start" in df.columns
    assert "end" in df.columns
    assert "strand" in df.columns


@task()
def load_enhancers(
    output_dir: str,
    tmpdir: str,
    candidate_enhancer_regions: str,
    chrom_sizes: str,
    features: dict,
    genes=None,
    force: bool = False,
    skip_rpkm_quantile: bool = False,
    celltype: str = None,
    tss_slop_for_class_assignment: int = 500,
    use_fast_count: bool = True,
    default_accessibility_feature: str = "",
    qnorm: str = None,
    class_override_file: str = None,
):
    makedirs(output_dir)
    makedirs(tmpdir)
    print("made dirs")
    local_candidate_enhancer_regions = join(
        tmpdir, basename(candidate_enhancer_regions)
    )
    download_file(candidate_enhancer_regions, tmpdir, overwrite_ok=True)
    enhancers = read_bed(local_candidate_enhancer_regions)
    enhancers = count_features_for_bed(
        enhancers,
        candidate_enhancer_regions,
        chrom_sizes,
        features,
        tmpdir,
        output_dir,
        "Enhancers",
        skip_rpkm_quantile,
        force,
        use_fast_count,
    )
    return process_enhancer_list(
        enhancers,
        celltype,
        genes,
        tss_slop_for_class_assignment,
        qnorm,
        default_accessibility_feature,
        tmpdir,
        output_dir,
    )


@task()
def process_enhancer_list(
    enhancers,
    celltype,
    genes,
    tss_slop_for_class_assignment,
    qnorm,
    default_accessibility_feature,
    tmpdir,
    output_dir,
):
    # celltype
    enhancers["chr"] = enhancers["chr"].astype("str")
    if celltype is not None:
        enhancers["celltype"] = celltype

    # Assign categories
    if genes is not None:
        print("Assigning classes to enhancers")
        enhancers = assign_enhancer_classes(
            enhancers, genes, tss_slop=tss_slop_for_class_assignment
        )

    # TO DO: Should qnorm each bam file separately (before averaging). Currently qnorm being performed on the average
    enhancers = run_qnorm(enhancers, qnorm)
    enhancers = compute_activity(enhancers, default_accessibility_feature)
    local_enhancer_list_out = join(tmpdir, "EnhancerList.txt")
    local_enhancer_list_bed_out = join(tmpdir, "EnhancerList.bed")
    enhancers.to_csv(
        local_enhancer_list_out, sep="\t", index=False, header=True, float_format="%.6f"
    )
    enhancers[["chr", "start", "end", "name"]].to_csv(
        local_enhancer_list_bed_out, sep="\t", index=False, header=False
    )

    # upload output files to remote output dir
    upload_file(local_enhancer_list_out, output_dir, overwrite_ok=True)
    upload_file(local_enhancer_list_bed_out, output_dir, overwrite_ok=True)
    return join(output_dir, "EnhancerList.txt"), join(output_dir, "EnhancerList.bed")


# Kristy's version
def assign_enhancer_classes(enhancers, genes, tss_slop=500):

    # build pyranges df
    tss_pyranges = df_to_pyranges(
        genes, start_col="tss", end_col="tss", start_slop=tss_slop, end_slop=tss_slop
    )
    gene_pyranges = df_to_pyranges(genes)

    def get_class_pyranges(
        enhancers, tss_pyranges=tss_pyranges, gene_pyranges=gene_pyranges
    ):
        """
        Takes in PyRanges objects : Enhancers, tss_pyranges, gene_pyranges
        Returns dataframe with  uid (representing enhancer) and symbol of the gene/promoter that is overlapped"""

        # genes
        genic_enh = enhancers.join(gene_pyranges, suffix="_genic")
        genic_enh = (
            genic_enh.df[["symbol", "uid"]]
            .groupby("uid", as_index=False)
            .aggregate(lambda x: ",".join(list(set(x))))
        )

        # promoters
        promoter_enh = enhancers.join(tss_pyranges, suffix="_promoter")
        promoter_enh = (
            promoter_enh.df[["symbol", "uid"]]
            .groupby("uid", as_index=False)
            .aggregate(lambda x: ",".join(list(set(x))))
        )

        return genic_enh, promoter_enh

    # label everything as intergenic
    enhancers["class"] = "intergenic"
    enhancers["uid"] = range(enhancers.shape[0])
    enh = df_to_pyranges(enhancers)

    genes, promoters = get_class_pyranges(enh)
    enhancers = enh.df.drop(["Chromosome", "Start", "End"], axis=1)
    enhancers.loc[enhancers["uid"].isin(genes.uid), "class"] = "genic"
    enhancers.loc[enhancers["uid"].isin(promoters.uid), "class"] = "promoter"

    enhancers["isPromoterElement"] = enhancers["class"] == "promoter"
    enhancers["isGenicElement"] = enhancers["class"] == "genic"
    enhancers["isIntergenicElement"] = enhancers["class"] == "intergenic"

    # Output stats
    print("Total enhancers: {}".format(len(enhancers)))
    print("         Promoters: {}".format(sum(enhancers["isPromoterElement"])))
    print("         Genic: {}".format(sum(enhancers["isGenicElement"])))
    print("         Intergenic: {}".format(sum(enhancers["isIntergenicElement"])))

    # Add promoter/genic symbol
    enhancers = enhancers.merge(
        promoters.rename(columns={"symbol": "promoterSymbol"}), on="uid", how="left"
    ).fillna(value={"promoterSymbol": ""})
    enhancers = enhancers.merge(
        genes.rename(columns={"symbol": "genicSymbol"}), on="uid", how="left"
    ).fillna(value={"genicSymbol": ""})
    enhancers.drop(["uid"], axis=1, inplace=True)

    # just to keep things consistent with original code
    enhancers["name"] = enhancers.apply(
        lambda e: "{}|{}:{}-{}".format(e["class"], e.chr, e.start, e.end), axis=1
    )
    return enhancers


@task()
def run_count_reads(
    target: str,
    output: str,
    output_dir: str,
    tmpdir: str,
    bed_file: str,
    chrom_sizes: str,
    use_fast_count: bool,
):
    if target.endswith(".bam"):
        return count_bam(
            bamfile=target,
            bed_file=bed_file,
            output=output,
            tmpdir=tmpdir,
            output_dir=output_dir,
            chrom_sizes=chrom_sizes,
            use_fast_count=use_fast_count,
        )
    elif target.endswith(".tagAlign.gz") or target.endswith(".tagAlign.bgz"):
        return count_tagalign(
            target, bed_file, output, output_dir, chrom_sizes, tmpdir=tmpdir
        )
    elif isBigWigFile(target):
        return count_bigwig(target, bed_file, output, output_dir)
    else:
        raise ValueError(
            "File {} name was not in .bam, .tagAlign.gz, .bw".format(target)
        )


@task()
def count_bam(
    bamfile: str,
    bed_file: str,
    output: str,
    tmpdir: str,
    output_dir: str,
    chrom_sizes: str,
    use_fast_count: bool = True,
    verbose: bool = True,
):
    local_bamfile = join(tmpdir, basename(bamfile))
    local_bed_file = join(tmpdir, basename(bed_file))
    local_chrom_sizes = join(tmpdir, basename(chrom_sizes))
    local_chrom_sizes_bed = ".".join([local_chrom_sizes, "bed"])
    chrom_sizes_bed = ".".join([chrom_sizes, "bed"])
    local_output = join(tmpdir, basename(output))
    if use_fast_count:
        temp_output = local_output + ".temp_sort_order"
        return script(
            f"""
        #!/bin/bash
        awk 'FNR==NR {{x2[$1] = $0; next}} $1 in x2 {{print x2[$1]}}' {local_chrom_sizes} <(samtools view -H {local_bamfile} | grep SQ | cut -f 2 | cut -c 4- )  > {temp_output};
        bedtools sort -faidx {temp_output} -i {local_bed_file} | bedtools coverage -g {temp_output} -counts  -a stdin -b {local_bamfile} | awk '{{print $1"\t"$2"\t"$3"\t"$NF}}'  | bedtools sort -faidx {local_chrom_sizes} -i stdin > {local_output}; rm {temp_output}
        """,
            inputs=[
                File(bamfile).stage(File(local_bamfile)),
                File(chrom_sizes).stage(File(local_chrom_sizes)),
                File(bed_file).stage(File(local_bed_file)),
            ],
            outputs={
                "file": File(join(output_dir, basename(output))).stage(
                    File(local_output)
                ),
                "path": join(output_dir, basename(output)),
            },
        )
    else:
        return script(
            f"""
        #!/bin/bash
        bedtools bamtobed -i {local_bamfile} | cut -f 1-3 | bedtools intersect -wa -a stdin -b {local_chrom_sizes_bed} | bedtools sort -i stdin -faidx {local_chrom_sizes} | bedtools coverage -g {local_chrom_sizes} -counts -sorted -a {local_bed_file} -b stdin | awk '{{print $1"\t"$2"\t"$3"\t"$NF}}' > {local_output}
        """,
            inputs=[
                File(bamfile).stage(File(local_bamfile)),
                File(chrom_sizes).stage(File(local_chrom_sizes)),
                File(chrom_sizes_bed).stage(File(local_chrom_sizes_bed)),
                File(bed_file).stage(File(local_bed_file)),
            ],
            outputs={
                "file": File(join(output_dir, basename(output))).stage(
                    File(local_output)
                ),
                "path": join(output_dir, basename(output)),
            },
        )


@task()
def count_tagalign(
    tagalign: str, bed_file: str, output: str, chrom_sizes: str, tmpdir: str
):
    return script(
        """
    #!/bin/bash
    tabix -B {tagalign} {bed_file} | cut -f1-3 |bedtools coverage -counts -b stdin -a {bed_file} | awk '{{print $1"\t"$2"\t" $3"\t"$NF}}'>{output}"""
    )


@task()
def count_bigwig(target: str, bed_file: str, output: str, tmpdir: str):
    from pyBigWig import open as open_bigwig

    bw = open_bigwig(target)
    bed = read_bed(bed_file)
    with open(output, "wb") as outfp:
        for chr, start, end, *rest in bed.itertuples(index=False, name=None):
            # if isinstance(name, np.float):
            #     name = ""
            try:
                val = (
                    bw.stats(chr, int(start), int(max(end, start + 1)), "mean")[0] or 0
                )
            except RuntimeError:
                print("Failed on", chr, start, end)
                raise
            val *= abs(end - start)  # convert to total coverage
            output = ("\t".join([chr, str(start), str(end), str(val)]) + "\n").encode(
                "ascii"
            )
            outfp.write(output)
    return output


def isBigWigFile(filename):
    return (
        filename.endswith(".bw")
        or filename.endswith(".bigWig")
        or filename.endswith(".bigwig")
    )


@task()
def count_features_for_bed(
    df,
    bed_file: str,
    chrom_sizes: str,
    features: str,
    tmpdir: str,
    output_dir: str,
    filebase: str,
    skip_rpkm_quantile: bool = False,
    force: bool = False,
    use_fast_count: bool = True,
):
    for feature, feature_bam_list in features.items():
        start_time = time.time()
        if isinstance(feature_bam_list, str):
            feature_bam_list = [feature_bam_list]

        for feature_bam in feature_bam_list:
            df = count_single_feature_for_bed(
                df,
                bed_file,
                chrom_sizes,
                feature_bam,
                feature,
                tmpdir,
                output_dir,
                filebase,
                skip_rpkm_quantile,
                force,
                use_fast_count,
            )

        df = average_features(
            df, feature.replace("feature_", ""), feature_bam_list, skip_rpkm_quantile
        )
        elapsed_time = time.time() - start_time
        print("Feature " + feature + " completed in " + str(elapsed_time))
    return df


@task()
def count_single_feature_for_bed(
    df,
    bed_file: str,
    chrom_sizes: str,
    feature_bam: str,
    feature: str,
    tmpdir: str,
    output_dir: str,
    filebase: str,
    skip_rpkm_quantile: bool,
    force: bool,
    use_fast_count: bool,
):
    orig_shape = df.shape[0]
    feature_name = feature + "." + basename(feature_bam)
    feature_outfile = join(
        output_dir, filebase + "." + feature_name + ".CountReads.bedgraph"
    )

    if (
        force
        or (not os.path.exists(feature_outfile))
        or (os.path.getsize(feature_outfile) == 0)
    ):
        print("Regenerating", feature_outfile)
        print("Counting coverage for {}".format(filebase + "." + feature_name))
        feature_outfile = run_count_reads(
            target=feature_bam,
            output=feature_outfile,
            output_dir=output_dir,
            tmpdir=tmpdir,
            bed_file=bed_file,
            chrom_sizes=chrom_sizes,
            use_fast_count=use_fast_count,
        )
    else:
        print(
            "Loading coverage from pre-calculated file for {}".format(
                filebase + "." + feature_name
            )
        )
    return process_feature_outfile(
        feature_outfile["path"],
        feature_bam,
        feature_name,
        tmpdir,
        df,
        orig_shape,
        skip_rpkm_quantile,
    )


@task()
def process_feature_outfile(
    feature_outfile,
    feature_bam,
    feature_name,
    tmpdir,
    df,
    orig_shape,
    skip_rpkm_quantile,
):
    if type(feature_outfile) == list:
        feature_outfile = feature_outfile[0]
    if type(feature_outfile) == File:
        feature_outfile = feature_outfile.path
    if feature_outfile.startswith("s3"):
        download_file(feature_outfile, tmpdir, overwrite_ok=True)
        feature_outfile = join(tmpdir, basename(feature_outfile))
    domain_counts = pd.read_csv(feature_outfile, sep="\t", header=None).rename(
        columns={0: "chr", 1: "start", 2: "end"}
    )  # read_bed(feature_outfile)
    score_column = domain_counts.columns[-1]

    total_counts = decode(count_total(feature_bam, tmpdir))
    domain_counts = domain_counts[["chr", "start", "end", score_column]]
    featurecount = feature_name + ".readCount"
    domain_counts.rename(columns={score_column: featurecount}, inplace=True)
    domain_counts["chr"] = domain_counts["chr"].astype("str")

    df = df.merge(domain_counts.drop_duplicates())
    # df = smart_merge(df, domain_counts.drop_duplicates())

    assert df.shape[0] == orig_shape, "Dimension mismatch"
    return annotate_feature_quantiles(
        df, total_counts, feature_name, skip_rpkm_quantile, featurecount
    )


@task()
def decode(counts):
    return int(counts.decode().strip())


@task()
def annotate_feature_quantiles(
    df, total_counts, feature_name, skip_rpkm_quantile, featurecount
):
    df[feature_name + ".RPM"] = 1e6 * df[featurecount] / float(total_counts)
    if not skip_rpkm_quantile:
        df[featurecount + ".quantile"] = df[featurecount].rank() / float(len(df))
        df[feature_name + ".RPM.quantile"] = df[feature_name + ".RPM"].rank() / float(
            len(df)
        )
        df[feature_name + ".RPKM"] = (
            1e3 * df[feature_name + ".RPM"] / (df.end - df.start).astype(float)
        )
        df[feature_name + ".RPKM.quantile"] = df[feature_name + ".RPKM"].rank() / float(
            len(df)
        )
    return df[~df.duplicated()]


@task()
def average_features(df, feature, feature_bam_list, skip_rpkm_quantile):
    feature_RPM_cols = [
        feature + "." + os.path.basename(feature_bam) + ".RPM"
        for feature_bam in feature_bam_list
    ]

    df[feature + ".RPM"] = df[feature_RPM_cols].mean(axis=1)

    if not skip_rpkm_quantile:
        feature_RPKM_cols = [
            feature + "." + os.path.basename(feature_bam) + ".RPKM"
            for feature_bam in feature_bam_list
        ]
        df[feature + ".RPM.quantile"] = df[feature + ".RPM"].rank() / float(len(df))
        df[feature + ".RPKM"] = df[feature_RPKM_cols].mean(axis=1)
        df[feature + ".RPKM.quantile"] = df[feature + ".RPKM"].rank() / float(len(df))

    return df


bed_extra_colnames = [
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
]
# JN: 9/13/19: Don't assume chromosomes start with 'chr'
# chromosomes = ['chr' + str(entry) for entry in list(range(1,23)) + ['M','X','Y']]   # should pass this in as an input file to specify chromosome order


@task()
def read_bed(
    filename,
    extra_colnames=bed_extra_colnames,
    chr=None,
    sort=False,
    skip_chr_sorting=True,
):
    skip = 1 if ("track" in open(filename, "r").readline()) else 0
    names = ["chr", "start", "end"] + extra_colnames
    result = pd.read_table(
        filename, names=names, header=None, skiprows=skip, comment="#"
    )
    result = result.dropna(axis=1, how="all")  # drop empty columns
    assert result.columns[0] == "chr"

    result["chr"] = pd.Categorical(result["chr"], ordered=True)
    if chr is not None:
        result = result[result.chr == chr]
    if not skip_chr_sorting:
        result.sort_values("chr", inplace=True)
    if sort:
        result.sort_values(["chr", "start", "end"], inplace=True)
    return result


@task()
def read_bedgraph(filename):
    read_bed(filename, extra_colnames=["score"], skip_chr_sorting=True)


@task()
def count_bam_mapped(bam_file, tmpdir):
    # Counts number of reads in a BAM file WITHOUT iterating.  Requires that the BAM is indexed
    if bam_file.startswith("s3"):
        local_bam_file = join(tmpdir, basename(bam_file))
        download_file(bam_file, tmpdir, overwrite_ok=True)
    else:
        local_bam_file = bam_file
    return script(
        f"""
    samtools index {local_bam_file}
    samtools idxstats {local_bam_file} | grep -v '*' |cut -f3 | paste -sd+ | bc 
    """
    )


@task()
def count_tagalign_total(tagalign, tmpdir):
    local_tagalign = join(tmpdir, basename(tagalign))
    return script(
        f"""
    zcat local_tagalign | grep -E 'chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY' | wc -l
    """,
        inputs=[File(tagalign).stage(local_tagalign)],
    )


@task()
def count_bigwig_total(bw_file, tmpdir):
    download_file(bw_file, tmpdir, overwrite_ok=True)
    bw_file = join(tmpdir, basename(bw_file))
    from pyBigWig import open as open_bigwig

    bw = open_bigwig(bw_file)
    result = sum(l * bw.stats(ch, 0, l, "mean")[0] for ch, l in bw.chroms().items())
    assert (
        abs(result) > 0
    )  ## BigWig could have negative values, e.g. the negative-strand GroCAP bigwigs
    return result


@task()
def count_total(infile, tmpdir):
    if infile.endswith(".tagAlign.gz") or infile.endswith(".tagAlign.bgz"):
        total_counts = count_tagalign_total(infile, tmpdir)
    elif infile.endswith(".bam"):
        total_counts = count_bam_mapped(infile, tmpdir)
    elif isBigWigFile(infile):
        total_counts = count_bigwig_total(infile, tmpdir)
    else:
        raise RuntimeError("Did not recognize file format of: " + infile)

    return total_counts


@task()
def parse_params_file(args):
    # Parse parameters file and return params dictionary
    params = {}

    params["default_accessibility_feature"] = determine_accessibility_feature(args)
    params["features"] = get_features(args)

    if args.expression_table:
        params["expression_table"] = args.expression_table.split(",")
    else:
        params["expression_table"] = ""

    return params


@task()
def get_features(h3k27ac, atac, dhs, supplementary_features):
    features = {}
    if h3k27ac:
        features["h3k27ac"] = h3k27ac.split(",")

    if atac:
        features["atac"] = atac.split(",")

    if dhs:
        features["dhs"] = dhs.split(",")

    if supplementary_features is not None:
        supp = pd.read_csv(supplementary_features, sep="\t")
        for idx, row in supp.iterrows():
            features[row["feature_name"]] = row["file"].split(",")
    print(str(features))
    return features


@task()
def determine_accessibility_feature(default_accessibility_feature, atac, dhs):
    if default_accessibility_feature is not None:
        return default_accessibility_feature
    elif (not atac) and (not dhs):
        raise RuntimeError(
            "Both dhs and atac have been provided. Must set one file to be the default accessibility feature!"
        )
    elif atac:
        return "atac"
    elif dhs:
        return "dhs"
    else:
        raise RuntimeError("At least one of atac or dhs must be provided!")


def compute_activity(df, access_col):
    if access_col == "dhs":
        if "h3k27ac.RPM" in df.columns:
            df["activity_base"] = np.sqrt(
                df["normalized_h3K27ac"] * df["normalized_dhs"]
            )
            df["activity_base_no_qnorm"] = np.sqrt(df["h3k27ac.RPM"] * df["dhs.RPM"])
        else:
            df["activity_base"] = df["normalized_dhs"]
            df["activity_base_no_qnorm"] = df["dhs.RPM"]
    elif access_col == "atac":
        if "h3k27ac.RPM" in df.columns:
            df["activity_base"] = np.sqrt(
                df["normalized_h3K27ac"] * df["normalized_atac"]
            )
            df["activity_base_no_qnorm"] = np.sqrt(df["h3k27ac.RPM"] * df["atac.RPM"])
        else:
            df["activity_base"] = df["normalized_atac"]
            df["activity_base_no_qnorm"] = df["atac.RPM"]
    else:
        raise RuntimeError("At least one of atac or dhs must be provided!")

    return df


def run_qnorm(df, qnorm, qnorm_method="rank", separate_promoters=True):
    # Quantile normalize epigenetic data to a reference
    #
    # Option to qnorm promoters and nonpromoters separately

    if qnorm is None:
        if "h3k27ac.RPM" in df.columns:
            df["normalized_h3K27ac"] = df["h3k27ac.RPM"]
        if "dhs.RPM" in df.columns:
            df["normalized_dhs"] = df["dhs.RPM"]
        if "atac.RPM" in df.columns:
            df["normalized_atac"] = df["atac.RPM"]
    else:
        qnorm = pd.read_csv(qnorm, sep="\t")
        nRegions = df.shape[0]
        col_dict = {
            "dhs.RPM": "normalized_dhs",
            "atac.RPM": "normalized_atac",
            "h3k27ac.RPM": "normalized_h3K27ac",
        }

        for col in set(df.columns & col_dict.keys()):
            # if there is no atac.RPM in the qnorm file, but there is atac.RPM in enhancers, then qnorm atac to dhs
            if col == "atac.RPM" and "atac.RPM" not in qnorm.columns:
                qnorm["atac.RPM"] = qnorm["dhs.RPM"]

            if not separate_promoters:
                qnorm = qnorm.loc[qnorm["enh_class" == "any"]]
                if qnorm_method == "rank":
                    interpfunc = interpolate.interp1d(
                        qnorm["rank"],
                        qnorm[col],
                        kind="linear",
                        fill_value="extrapolate",
                    )
                    df[col_dict[col]] = interpfunc(
                        (1 - df[col + ".quantile"]) * nRegions
                    ).clip(0)
                elif qnorm_method == "quantile":
                    interpfunc = interpolate.interp1d(
                        qnorm["quantile"],
                        qnorm[col],
                        kind="linear",
                        fill_value="extrapolate",
                    )
                    df[col_dict[col]] = interpfunc(df[col + ".quantile"]).clip(0)
            else:
                for enh_class in ["promoter", "nonpromoter"]:
                    this_qnorm = qnorm.loc[qnorm["enh_class"] == enh_class]

                    # Need to recompute quantiles within each class
                    if enh_class == "promoter":
                        this_idx = df.index[
                            np.logical_or(
                                df["class"] == "tss", df["class"] == "promoter"
                            )
                        ]
                    else:
                        this_idx = df.index[
                            np.logical_and(
                                df["class"] != "tss", df["class"] != "promoter"
                            )
                        ]
                    df.loc[this_idx, col + enh_class + ".quantile"] = df.loc[
                        this_idx, col
                    ].rank() / len(this_idx)

                    if qnorm_method == "rank":
                        interpfunc = interpolate.interp1d(
                            this_qnorm["rank"],
                            this_qnorm[col],
                            kind="linear",
                            fill_value="extrapolate",
                        )
                        df.loc[this_idx, col_dict[col]] = interpfunc(
                            (1 - df.loc[this_idx, col + enh_class + ".quantile"])
                            * len(this_idx)
                        ).clip(0)
                    elif qnorm_method == "quantile":
                        interpfunc = interpolate.interp1d(
                            this_qnorm["quantile"],
                            this_qnorm[col],
                            kind="linear",
                            fill_value="extrapolate",
                        )
                        df.loc[this_idx, col_dict[col]] = interpfunc(
                            df.loc[this_idx, col + enh_class + ".quantile"]
                        ).clip(0)

    return df
