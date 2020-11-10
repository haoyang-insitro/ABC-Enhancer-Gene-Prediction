import pdb
import sys
import time
import traceback
from os.path import basename, join

import numpy as np
import pandas as pd
from predictor import *
from redun import Dir, File, script, task
from tools import *

from insitro_core.utils.cloud.bucket_utils import download_file
from insitro_core.utils.storage import *


@task()
def predict(
    output_dir: str,
    tmpdir: str,
    enhancers: str,
    genes: str,
    celltype: str,
    hic_resolution: int,
    window: int = 5e6,
    score_column: str = "ABC.Score",
    threshold: float = 0.022,
    hicdir: str = None,
    tss_hic_contribution: float = 100,
    hic_pseudocount_distance: int = 1e6,
    hic_type: str = "juicebox",
    hic_is_doubly_stochastic: bool = False,
    scale_hic_using_powerlaw: bool = False,
    hic_gamma: float = 0.87,
    hic_gamma_reference: float = 0.87,
    run_all_genes: bool = False,
    expression_cutoff: float = 1,
    promoter_activity_quantile_cutoff: float = 0.4,
    make_all_putative: bool = False,
    use_hdf5: bool = False,
    tss_slop: float = 500,
    chromosomes: str = "all",
    include_chrY: bool = True,
):
    """
    #Default Params:
    enhancers: Candidate enhancer regions. Formatted as the EnhancerList.txt file produced by run.neighborhoods.py
    genes: Genes to make predictions for. Formatted as the GeneList.txt file produced by run.neighborhoods.py
    celltype: Name of cell type.
    window: Make predictions for all candidate elements within this distance of the gene's TSS.
    score_column: Column name of score to use for thresholding.
    threshold: Threshold on ABC Score (--score_column) to call a predicted positive

    #HiC params:
    hicdir: HiC directory.
    hic_resolution: HiC resolution.
    tss_hic_contribution: Weighting of diagonal bin of hic matrix as a percentage of the maximum of its neighboring bins
    hic_pseudocount_distance:A pseudocount is added equal to the powerlaw fit at this distance
    hic_type:format of hic files; one of juicebox, bedpe
    hich_is_doubly_stochastic: If hic matrix is already DS, can skip this step

    #PowerLaw params:
    scale_hic_using_powerlaw: Scale Hi-C values using powerlaw relationship
    hic_gamma: Powerlaw exponent of hic data. Must be positive
    run_all_genes: Do not check for gene expression, make predictions for all genes
    expression_cutoff: Make predictions for genes with expression higher than this value
    promoter_activity_quantile_cutoff: Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data
    make_all_putative: Make big file with concatenation of all genes file
    use_hdf5: Write AllPutative file in hdf5 format instead of tab-delimited
    tss_slop: Distance from tss to search for self-promoters
    chromosomes: Chromosomes to make predictions for. Defaults to intersection of all chromosomes in --genes and --enhancers
    inlcude_chrY: Make predictions on Y chromosome
    """
    makedirs(output_dir)
    makedirs(tmpdir)

    # write the parameters
    write_params(
        {
            "output_dir": output_dir,
            "tmpdir": tmpdir,
            "enhancers": enhancers,
            "genes": genes,
            "celltype": celltype,
            "hic_resolution": hic_resolution,
            "window": window,
            "score_column": score_column,
            "threshold": threshold,
            "hicdir": hicdir,
            "tss_hic_contribution": tss_hic_contribution,
            "hic_pseudocount_distance": hic_pseudocount_distance,
            "hic_type": hic_type,
            "hic_is_doubly_stochastic": hic_is_doubly_stochastic,
            "scale_hic_using_powerlaw": scale_hic_using_powerlaw,
            "hic_gamma": hic_gamma,
            "hic_gamma_reference": hic_gamma_reference,
            "run_all_genes": run_all_genes,
            "expression_cutoff": expression_cutoff,
            "promoter_activity_quantile_cutoff": promoter_activity_quantile_cutoff,
            "make_all_putative": make_all_putative,
            "use_hdf5": use_hdf5,
            "tss_slop": tss_slop,
            "chromosomes": chromosomes,
            "include_chrY": include_chrY,
        },
        output_dir=output_dir,
        tmpdir=tmpdir,
        output_fname="parameters.predict.txt",
    )

    if type(genes) == str:
        print("reading genes")
        genes = pd.read_csv(genes, sep="\t")
    genes = determine_expressed_genes(
        genes, expression_cutoff, promoter_activity_quantile_cutoff
    )
    genes = genes.loc[
        :,
        [
            "chr",
            "symbol",
            "tss",
            "Expression",
            "PromoterActivityQuantile",
            "isExpressed",
        ],
    ]
    genes.columns = [
        "chr",
        "TargetGene",
        "TargetGeneTSS",
        "TargetGeneExpression",
        "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed",
    ]

    print("reading enhancers")
    enhancers_full = pd.read_csv(enhancers, sep="\t")
    # TO DO
    # Think about which columns to include
    enhancers = enhancers_full.loc[
        :, ["chr", "start", "end", "name", "class", "activity_base"]
    ]

    # Initialize Prediction files
    pred_file_full = join(tmpdir, "EnhancerPredictionsFull.txt")
    pred_file_slim = join(tmpdir, "EnhancerPredictions.txt")
    pred_file_bedpe = join(tmpdir, "EnhancerPredictions.bedpe")
    all_pred_file_expressed = join(tmpdir, "EnhancerPredictionsAllPutative.txt.gz")
    all_pred_file_nonexpressed = join(
        tmpdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz"
    )
    all_putative_list = []

    # Make predictions
    if chromosomes == "all":
        chromosomes = set(genes["chr"]).intersection(set(enhancers["chr"]))
        if not include_chrY:
            chromosomes.discard("chrY")
    else:
        chromosomes = chromosomes.split(",")

    for chromosome in chromosomes:
        print("Making predictions for chromosome: {}".format(chromosome))
        t = time.time()

        this_enh = enhancers.loc[enhancers["chr"] == chromosome, :].copy()
        this_genes = genes.loc[genes["chr"] == chromosome, :].copy()

        this_chr = make_predictions(
            tmpdir,
            chromosome,
            this_enh,
            this_genes,
            window,
            tss_slop,
            hic_type,
            hicdir,
            hic_gamma,
            hic_gamma_reference,
            hic_resolution,
            tss_hic_contribution,
            scale_hic_using_powerlaw,
            hic_pseudocount_distance,
        )
        all_putative_list.append(this_chr)

        print(
            "Completed chromosome: {}. Elapsed time: {} \n".format(
                chromosome, time.time() - t
            )
        )

    # Subset predictions
    print("Writing output files...")
    all_putative = pd.concat(all_putative_list)
    all_putative["CellType"] = celltype
    slim_cols = [
        "chr",
        "start",
        "end",
        "name",
        "TargetGene",
        "TargetGeneTSS",
        "CellType",
        score_column,
    ]
    if run_all_genes:
        all_positive = all_putative.iloc[
            np.logical_and.reduce(
                (
                    all_putative[score_column] > threshold,
                    ~(all_putative["class"] == "promoter"),
                )
            ),
            :,
        ]
    else:
        all_positive = all_putative.iloc[
            np.logical_and.reduce(
                (
                    all_putative.TargetGeneIsExpressed,
                    all_putative[score_column] > threshold,
                    ~(all_putative["class"] == "promoter"),
                )
            ),
            :,
        ]

    all_positive.to_csv(
        pred_file_full, sep="\t", index=False, header=True, float_format="%.6f"
    )
    all_positive[slim_cols].to_csv(
        pred_file_slim, sep="\t", index=False, header=True, float_format="%.6f"
    )

    make_gene_prediction_stats(all_putative, score_column, threshold, tmpdir)
    write_connections_bedpe_format(all_positive, pred_file_bedpe, score_column)

    if make_all_putative:
        if not use_hdf5:
            all_putative.loc[all_putative.TargetGeneIsExpressed, :].to_csv(
                all_pred_file_expressed,
                sep="\t",
                index=False,
                header=True,
                compression="gzip",
                float_format="%.6f",
                na_rep="NaN",
            )
            all_putative.loc[~all_putative.TargetGeneIsExpressed, :].to_csv(
                all_pred_file_nonexpressed,
                sep="\t",
                index=False,
                header=True,
                compression="gzip",
                float_format="%.6f",
                na_rep="NaN",
            )
        else:
            all_pred_file_expressed = join(tmpdir, "EnhancerPredictionsAllPutative.h5")
            all_pred_file_nonexpressed = join(
                tmpdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.h5"
            )
            all_putative.loc[all_putative.TargetGeneIsExpressed, :].to_hdf(
                all_pred_file_expressed, key="predictions", complevel=9, mode="w"
            )
            all_putative.loc[~all_putative.TargetGeneIsExpressed, :].to_hdf(
                all_pred_file_nonexpressed, key="predictions", complevel=9, mode="w"
            )

    # copy local files to remote
    upload_file(pred_file_full, output_dir, overwrite_ok=True)
    upload_file(pred_file_slim, output_dir, overwrite_ok=True)
    upload_file(pred_file_bedpe, output_dir, overwrite_ok=True)
    upload_file(all_pred_file_expressed, output_dir, overwrite_ok=True)
    upload_file(all_pred_file_nonexpressed, output_dir, overwrite_ok=True)
    print("Done.")
    return [
        pred_file_full,
        pred_file_slim,
        pred_file_bedpe,
        all_pred_file_expressed,
        all_pred_file_nonexpressed,
    ]
