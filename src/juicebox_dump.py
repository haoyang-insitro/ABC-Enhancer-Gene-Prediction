from os.path import basename, join

import redun
from redun import Dir, script, task

from insitro_core.utils.cloud.bucket_utils import download_file
from insitro_core.utils.storage import *


@task()
def download_observed_matrix(
    juicebox: str,
    hic_file: str,
    chromosome: str,
    output_dir: str,
    resolution: int,
    tmpdir: str,
):
    return script(
        f"""
    #!/bin/bash 
    {juicebox} dump observed KR {hic_file} {chromosome} {chromosome} BP {resolution} {tmpdir}/chr{chromosome}.KRobserved
    gzip -f {tmpdir}/chr{chromosome}.KRobserved
    {juicebox} dump norm KR {hic_file} {chromosome} BP {resolution} {tmpdir}/chr{chromosome}.KRnorm
    gzip -f {tmpdir}/chr{chromosome}.KRnorm
    """,
        outputs=[Dir(output_dir).stage(Dir(tmpdir))],
    )


@task()
def download_raw(
    juicebox: str,
    hic_file: str,
    chromosome: str,
    output_dir: str,
    resolution: int,
    tmpdir: str,
):
    return script(
        f"""
    #!/bin/bash 
    {juicebox} dump observed NONE {hic_file} {chromosome} {chromosome} BP {resolution} {tmpdir}/chr{chromosome}.RAWobserved
    gzip -f {tmpdir}/chr{chromosome}.RAWobserved
    """,
        outputs=[Dir(output_dir).stage(Dir(tmpdir))],
    )


@task()
def download_hic(
    output_dir: str,
    tmpdir: str,
    hic_file: str,
    juicebox: str = "",
    resolution: int = 5000,
    include_raw: bool = False,
    chromosomes: str = "all",
    skip_gzip: bool = False,
):
    """
    #Params
    hic_file: Path or url to .hic file
    juicebox: path to juicebox executable or java command invoking juicer_tools.jar. eg: 'java -jar juicer_tools.jar'
    resolution: Resolution of HiC to download. In units of bp.
    include_raw: Download raw matrix in addition to KR
    chromosomes: comma delimited list of chromosomes to download
    skip_gzip: dont gzip hic files
    """
    makedirs(output_dir)
    makedirs(tmpdir)

    if chromosomes == "all":
        chromosomes = list(range(1, 23)) + ["X"]
    else:
        chromosomes = chromosomes.split(",")

    for chromosome in chromosomes:
        chrom_tmpdir = join(tmpdir, "chr" + chromosome)
        chrom_output_dir = join(output_dir, "chr" + chromosome)
        makedirs(chrom_tmpdir)
        makedirs(chrom_output_dir)

        result = download_observed_matrix(
            juicebox, hic_file, chromosome, chrom_output_dir, resolution, chrom_tmpdir
        )
        if include_raw:
            result_raw = download_raw(
                juicebox,
                hic_file,
                chromosome,
                chrom_output_dir,
                resolution,
                chrom_tmpdir,
            )
            return result, result_raw
        else:
            return result
