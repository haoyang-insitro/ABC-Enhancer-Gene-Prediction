import glob
import sys
import traceback
from os.path import basename, join

import numpy as np
import pandas
from hic import *
from redun import Dir, task
from scipy import stats

from insitro_core.utils.cloud.bucket_utils import *
from insitro_core.utils.storage import *


@task()
def compute_powerlaw_fit_from_hic(
    hicdir: str,
    output_dir: str,
    tmpdir: str,
    hic_type: str = "juicebox",
    resolution: int = 5000,
    minwindow: int = 5000,
    maxwindow: int = 1e6,
    chromosomes: str = "all",
):
    """
    #Params
    hicdir:Directory containing observed HiC KR normalized matrices. File naming and structure should be: hicDir/chr*/chr*.KRobserved
    hic_type: format of hic files
    resolution:For Juicebox: resolution of hic dataset (in bp). For bedpe: distances will be binned to this resolution for powerlaw fit
    minwindow:Minimum distance between bins to include in powerlaw fit (bp). Recommended to be at least >= resolution to avoid using the diagonal of the HiC Matrix
    maxwindow:Maximum distance between bins to include in powerlaw fit (bp)
    chromosomes:Comma delimited list of chromosomes to use for fit. Defualts to chr[1..22],chrX
    """
    makedirs(output_dir)
    HiC = load_hic_for_powerlaw(
        tmpdir, chromosomes, hic_type, hicdir, resolution, minwindow, maxwindow
    )

    # Run
    slope, intercept, hic_mean_var = do_powerlaw_fit(HiC)

    # Store outputs in outputdir
    res = pandas.DataFrame(
        {
            "resolution": [resolution],
            "maxWindow": [maxwindow],
            "minWindow": [minwindow],
            "pl_gamma": [slope],
            "pl_scale": [intercept],
        }
    )
    res.to_csv(join(tmpdir, "hic.powerlaw.txt"), sep="\t", index=False, header=True)
    hic_mean_var.to_csv(
        join(tmpdir, "hic.mean_var.txt"), sep="\t", index=True, header=True
    )
    upload_file(join(tmpdir, "hic.powerlaw.txt"), output_dir)
    upload_file(join(tmpdir, "hic.mean_var.txt"), output_dir)


def load_hic_for_powerlaw(
    tmpdir, chromosomes, hic_type, hicdir, resolution, minwindow, maxwindow
):
    if chromosomes == "all":
        chromosomes = ["chr" + str(x) for x in list(range(1, 23))] + ["chrX"]
    else:
        chromosomes = chromosomes.split(",")

    all_data_list = []
    for chrom in chromosomes:
        try:
            if hic_type == "juicebox":
                hic_file, hic_norm_file, hic_is_vc = get_hic_file(
                    chrom, hicdir, allow_vc=False
                )
                print("Working on {}".format(hic_file))
                this_data = load_hic(
                    tmpdir=tmpdir,
                    hic_file=hic_file,
                    hic_norm_file=hic_norm_file,
                    hic_is_vc=hic_is_vc,
                    hic_type="juicebox",
                    hic_resolution=resolution,
                    tss_hic_contribution=100,
                    window=maxwindow,
                    min_window=minwindow,
                    gamma=np.nan,
                    interpolate_nan=False,
                )
                this_data["dist_for_fit"] = (
                    abs(this_data["bin1"] - this_data["bin2"]) * resolution
                )
                all_data_list.append(this_data)
            elif hic_type == "bedpe":
                hic_file = get_hic_file(chrom, hicdir, hic_type="bedpe")
                print("Working on {}".format(hic_file))
                this_data = load_hic(
                    hic_file=hic_file,
                    tmpdir=tmpdir,
                    hic_type="bedpe",
                    hic_norm_file=None,
                    hic_is_vc=None,
                    hic_resolution=None,
                    tss_hic_contribution=None,
                    window=None,
                    min_window=None,
                    gamma=None,
                )

                # Compute distance in bins as with juicebox data.
                # This is needed to in order to maintain consistency, but is probably slightly less accurate.
                # Binning also reduces noise level.
                rawdist = abs(
                    (this_data["x2"] + this_data["x1"]) / 2
                    - (this_data["y2"] + this_data["y1"]) / 2
                )
                this_data["dist_for_fit"] = (rawdist // resolution) * resolution
                this_data = this_data.loc[
                    np.logical_and(
                        this_data["dist_for_fit"] >= minwindow,
                        this_data["dist_for_fit"] <= maxwindow,
                    )
                ]
                all_data_list.append(this_data)
            else:
                error("invalid --hic_type")

        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)

    all_data = pd.concat(all_data_list)
    return all_data


def do_powerlaw_fit(HiC):
    print("Running regression")

    # TO DO:
    # Print out mean/var plot of powerlaw relationship
    HiC_summary = HiC.groupby("dist_for_fit").agg({"hic_contact": "sum"})
    HiC_summary["hic_contact"] = (
        HiC_summary.hic_contact / HiC_summary.hic_contact.sum()
    )  # technically this normalization should be over the entire genome (not just to maxWindow). Will only affect intercept though..
    res = stats.linregress(
        np.log(HiC_summary.index), np.log(HiC_summary["hic_contact"])
    )

    hic_mean_var = HiC.groupby("dist_for_fit").agg({"hic_contact": ["mean", "var"]})
    hic_mean_var.columns = ["mean", "var"]

    return res.slope, res.intercept, hic_mean_var
