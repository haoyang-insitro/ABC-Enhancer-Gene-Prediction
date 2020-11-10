import os
import pdb
import sys
import time

import numpy as np
import pandas as pd
import pyranges as pr
from hic import *
from tools import *


def make_predictions(
    tmpdir,
    chromosome,
    enhancers,
    genes,
    window,
    tss_slop,
    hic_type,
    HiCdir,
    hic_gamma,
    hic_gamma_reference,
    hic_resolution,
    tss_hic_contribution,
    scale_hic_using_powerlaw,
    hic_pseudocount_distance,
):
    pred = make_pred_table(chromosome, enhancers, genes, window)
    pred = annotate_predictions(pred, tss_slop)
    pred = add_powerlaw_to_predictions(pred, hic_gamma, hic_gamma_reference)

    # if Hi-C directory is not provided, only powerlaw model will be computed
    if HiCdir:
        hic_file, hic_norm_file, hic_is_vc = get_hic_file(
            chromosome, HiCdir, hic_type=hic_type
        )
        pred = add_hic_to_enh_gene_table(
            tmpdir,
            enhancers,
            genes,
            pred,
            hic_file,
            hic_norm_file,
            hic_is_vc,
            chromosome,
            hic_type,
            hic_resolution,
            tss_hic_contribution,
            window,
            hic_gamma,
            scale_hic_using_powerlaw,
            hic_pseudocount_distance,
        )
        pred = compute_score(
            pred, [pred["activity_base"], pred["hic_contact_pl_scaled_adj"]], "ABC"
        )

    pred = compute_score(
        pred, [pred["activity_base"], pred["powerlaw_contact_reference"]], "powerlaw"
    )

    return pred


def make_pred_table(chromosome, enh, genes, window):
    print("Making putative predictions table...")
    t = time.time()
    enh["enh_midpoint"] = (enh["start"] + enh["end"]) / 2
    enh["enh_idx"] = enh.index
    genes["gene_idx"] = genes.index
    enh_pr = df_to_pyranges(enh)
    genes_pr = df_to_pyranges(
        genes,
        start_col="TargetGeneTSS",
        end_col="TargetGeneTSS",
        start_slop=window,
        end_slop=window,
    )

    pred = enh_pr.join(genes_pr).df.drop(
        ["Start_b", "End_b", "chr_b", "Chromosome", "Start", "End"], axis=1
    )
    pred["distance"] = abs(pred["enh_midpoint"] - pred["TargetGeneTSS"])
    pred = pred.loc[pred["distance"] < window, :]  # for backwards compatability
    return pred


def add_hic_to_enh_gene_table(
    tmpdir,
    enh,
    genes,
    pred,
    hic_file,
    hic_norm_file,
    hic_is_vc,
    chromosome,
    hic_type,
    hic_resolution,
    tss_hic_contribution,
    window,
    hic_gamma,
    scale_hic_using_powerlaw,
    hic_pseudocount_distance,
):
    print("Begin HiC")
    HiC = load_hic(
        tmpdir=tmpdir,
        hic_file=hic_file,
        hic_norm_file=hic_norm_file,
        hic_is_vc=hic_is_vc,
        hic_type=hic_type,
        hic_resolution=hic_resolution,
        tss_hic_contribution=tss_hic_contribution,
        window=window,
        min_window=0,
        gamma=hic_gamma,
    )

    # Add hic to pred table
    # At this point we have a table where each row is an enhancer/gene pair.
    # We need to add the corresponding HiC matrix entry.
    # If the HiC is provided in juicebox format (ie constant resolution), then we can just merge using the indices
    # But more generally we do not want to assume constant resolution. In this case hic should be provided in bedpe format

    t = time.time()
    if hic_type == "bedpe":
        # Use pyranges to compute overlaps between enhancers/genes and hic bedpe table
        # Consider each range of the hic matrix separately - and merge each range into both enhancers and genes.
        # Then remerge on hic index

        HiC["hic_idx"] = HiC.index
        hic1 = df_to_pyranges(HiC, start_col="x1", end_col="x2", chr_col="chr1")
        hic2 = df_to_pyranges(HiC, start_col="y1", end_col="y2", chr_col="chr2")

        # Overlap in one direction
        enh_hic1 = (
            df_to_pyranges(
                enh, start_col="enh_midpoint", end_col="enh_midpoint", end_slop=1
            )
            .join(hic1)
            .df
        )
        genes_hic2 = (
            df_to_pyranges(
                genes, start_col="TargetGeneTSS", end_col="TargetGeneTSS", end_slop=1
            )
            .join(hic2)
            .df
        )
        ovl12 = enh_hic1[["enh_idx", "hic_idx", "hic_contact"]].merge(
            genes_hic2[["gene_idx", "hic_idx"]], on="hic_idx"
        )

        # Overlap in the other direction
        enh_hic2 = (
            df_to_pyranges(
                enh, start_col="enh_midpoint", end_col="enh_midpoint", end_slop=1
            )
            .join(hic2)
            .df
        )
        genes_hic1 = (
            df_to_pyranges(
                genes, start_col="TargetGeneTSS", end_col="TargetGeneTSS", end_slop=1
            )
            .join(hic1)
            .df
        )
        ovl21 = enh_hic2[["enh_idx", "hic_idx", "hic_contact"]].merge(
            genes_hic1[["gene_idx", "hic_idx"]], on=["hic_idx"]
        )

        # Concatenate both directions and merge into preditions
        ovl = pd.concat([ovl12, ovl21]).drop_duplicates()
        pred = pred.merge(ovl, on=["enh_idx", "gene_idx"], how="left")
        pred.fillna(value={"hic_contact": 0}, inplace=True)
    elif hic_type == "juicebox":
        # Merge directly using indices
        # Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
        # Index into sparse matrix
        # pred['hic_contact'] = [HiC[i,j] for (i,j) in pred[['enh_bin','tss_bin']].values.tolist()]

        pred["enh_bin"] = np.floor(pred["enh_midpoint"] / hic_resolution).astype(int)
        pred["tss_bin"] = np.floor(pred["TargetGeneTSS"] / hic_resolution).astype(int)
        if not hic_is_vc:
            # in this case the matrix is upper triangular.
            #
            pred["bin1"] = np.amin(pred[["enh_bin", "tss_bin"]], axis=1)
            pred["bin2"] = np.amax(pred[["enh_bin", "tss_bin"]], axis=1)
            pred = pred.merge(HiC, how="left", on=["bin1", "bin2"])
            pred.fillna(value={"hic_contact": 0}, inplace=True)
        else:
            # The matrix is not triangular, its full
            # For VC assume genes correspond to rows and columns to enhancers
            pred = pred.merge(
                HiC,
                how="left",
                left_on=["tss_bin", "enh_bin"],
                right_on=["bin1", "bin2"],
            )

        pred.fillna(value={"hic_contact": 0}, inplace=True)

        # QC juicebox HiC
        pred = qc_hic(pred)

    pred.drop(
        [
            "x1",
            "x2",
            "y1",
            "y2",
            "bin1",
            "bin2",
            "enh_idx",
            "gene_idx",
            "hic_idx",
            "enh_midpoint",
            "tss_bin",
            "enh_bin",
        ],
        inplace=True,
        axis=1,
        errors="ignore",
    )

    print("HiC added to predictions table. Elapsed time: {}".format(time.time() - t))

    # Add powerlaw scaling
    pred = scale_hic_with_powerlaw(pred, scale_hic_using_powerlaw)

    # Add pseudocount
    pred = add_hic_pseudocount(pred, hic_gamma, hic_pseudocount_distance)

    print("HiC Complete")
    # print('Elapsed time: {}'.format(time.time() - t))

    return pred


def scale_hic_with_powerlaw(pred, scale_hic_using_powerlaw):
    # Scale hic values to reference powerlaw
    if not scale_hic_using_powerlaw:
        pred["hic_contact_pl_scaled"] = pred["hic_contact"]
    else:
        pred["hic_contact_pl_scaled"] = pred["hic_contact"] * (
            pred["powerlaw_contact_reference"] / pred["powerlaw_contact"]
        )

    return pred


def add_powerlaw_to_predictions(pred, hic_gamma, hic_gamma_reference):
    pred["powerlaw_contact"] = get_powerlaw_at_distance(
        pred["distance"].values, hic_gamma
    )
    pred["powerlaw_contact_reference"] = get_powerlaw_at_distance(
        pred["distance"].values, hic_gamma_reference
    )

    return pred


def add_hic_pseudocount(pred, hic_gamma, hic_pseudocount_distance):
    # Add a pseudocount based on the powerlaw expected count at a given distance

    powerlaw_fit = get_powerlaw_at_distance(pred["distance"].values, hic_gamma)
    powerlaw_fit_at_ref = get_powerlaw_at_distance(hic_pseudocount_distance, hic_gamma)

    pseudocount = np.amin(
        pd.DataFrame({"a": powerlaw_fit, "b": powerlaw_fit_at_ref}), axis=1
    )
    pred["hic_pseudocount"] = pseudocount
    pred["hic_contact_pl_scaled_adj"] = pred["hic_contact_pl_scaled"] + pseudocount

    return pred


def qc_hic(pred, threshold=0.01):
    # Genes with insufficient hic coverage should get nan'd

    summ = (
        pred.loc[pred["isSelfPromoter"], :]
        .groupby(["TargetGene"])
        .agg({"hic_contact": "sum"})
    )
    bad_genes = summ.loc[summ["hic_contact"] < threshold, :].index

    pred.loc[pred["TargetGene"].isin(bad_genes), "hic_contact"] = np.nan

    return pred


def compute_score(enhancers, product_terms, prefix):

    scores = np.column_stack(product_terms).prod(axis=1)

    enhancers[prefix + ".Score.Numerator"] = scores
    enhancers[prefix + ".Score"] = enhancers[
        prefix + ".Score.Numerator"
    ] / enhancers.groupby("TargetGene")[prefix + ".Score.Numerator"].transform("sum")

    return enhancers


def annotate_predictions(pred, tss_slop=500):
    # TO DO: Add is self genic
    pred["isSelfPromoter"] = np.logical_and.reduce(
        (
            pred["class"] == "promoter",
            pred.start - tss_slop < pred.TargetGeneTSS,
            pred.end + tss_slop > pred.TargetGeneTSS,
        )
    )

    return pred


def make_gene_prediction_stats(pred, score_column, threshold, tmpdir):
    summ1 = pred.groupby(["chr", "TargetGene", "TargetGeneTSS"]).agg(
        {
            "TargetGeneIsExpressed": lambda x: set(x).pop(),
            score_column: lambda x: all(np.isnan(x)),
            "name": "count",
        }
    )
    summ1.columns = ["geneIsExpressed", "geneFailed", "nEnhancersConsidered"]

    summ2 = (
        pred.loc[pred["class"] != "promoter", :]
        .groupby(["chr", "TargetGene", "TargetGeneTSS"])
        .agg({score_column: lambda x: sum(x > threshold)})
    )
    summ2.columns = ["nDistalEnhancersPredicted"]
    if (summ1.shape[0] > 0) and (summ2.shape[0]) > 0:
        summ1 = summ1.merge(summ2, left_index=True, right_index=True)
        summ1.to_csv(
            os.path.join(tmpdir, "GenePredictionStats.txt"), sep="\t", index=True
        )
    # else: write empty outfile
    elif summ1.shape[0] == 0:
        summ1.to_csv(
            os.path.join(tmpdir, "GenePredictionStats.txt"), sep="\t", index=True
        )
    else:
        summ2.to_csv(
            os.path.join(tmpdir, "GenePredictionStats.txt"), sep="\t", index=True
        )
