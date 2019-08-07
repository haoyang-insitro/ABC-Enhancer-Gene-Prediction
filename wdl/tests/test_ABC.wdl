version 1.0

import "../wdl/abc-pipeline.wdl" as target
import "ValidateABC.wdl" as checker

workflow TestAbcPR {
   input {
        File dnaseqbam
        File dnaseqbam_index
        File chrom_sizes
        File? regions_blacklist
        File? regions_whitelist
        File genes_bed
        File h3k27ac_bam
        File h3k27ac_bam_index
        File? expression_table
        File? ubiq_genes
        File HiCdirTar
        String cellType = "defCellType"
        Float threshold = 0.022
        String candidateRegions_targetMD5
   }

   call target.ABCpipeline as target {
     input:
         dnaseqbam = dnaseqbam,
         dnaseqbam_index = dnaseqbam_index,
         dhs_bam = dnaseqbam,
         dhs_bam_index = dnaseqbam_index,
         chrom_sizes = chrom_sizes,
         regions_blacklist = regions_blacklist,
         regions_whitelist = regions_whitelist,
         genes_bed = genes_bed,
         h3k27ac_bam = h3k27ac_bam,
         h3k27ac_bam_index = h3k27ac_bam_index,
         expression_table = expression_table,
         ubiq_genes = ubiq_genes,
         HiCdirTar = HiCdirTar,
         cellType = cellType,
         threshold = threshold
   }

   call checker.ValidateABC as checker {
      input:
         candidateRegions_targetMD5 = candidateRegions_targetMD5,
         candidateRegions = target.candidateRegions
   }
}
