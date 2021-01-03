#' Model-based analysis of alternative polyadenylation using 3’ end-linked reads
#'
#' @param gtf A character specifying the full path of the GTF file (reference genome);
#' @param pas_annotation A list containing the pas annotation. MAAPER provides processed annotation
#' information from PolyA_DB v3 on its Github page.
#' @param output_dir A character specifying the full path of the output directory,
#' which is used to store all intermdediate and final outputs.
#' @param bam_c1 A character vector specifying the full paths to the bam files for
#' condition 1 (control). The length of the vector equals the number of samples.
#' @param bam_c2 A character vector specifying the full paths to the bam files for
#' condition 2 (experiment). The length of the vector equals the number of samples.
#' @param read_len An integer specifying the read length.
#' @param ncores An integer specifying the number of cores used in parallel computation.
#' @param num_pas_thre An integer specifying the threhold on PAS's read number. Defaults to 25.
#' @param frac_pas_thre A numeric specifying the threshold on PAS's fraction. Defaults to 0.05.
#' @param dist_thre An integer specifying the threshold on fragment length. Defaults to 600.
#' @param num_thre An integer specifying the threhold on gene's read number. Defaults to 50.
#' @param run "all" (default) or "skip-train". For test and debug only.
#' @param subset A character vector specifying genes' Ensembl IDs if only a subset of genes need to be analyzed.
#' Check the \code{pas_annotation} files for ID formats.
#' @param region "all" (default). For test and debug only.
#' @param gtf_rds NULL (default). For test and debug only.
#' @param verbose FALSE (default). For test and debug only.
#' @return \code{mapper} saves two text files, gene.txt and pas.txt, to \code{out_dir}.
#' pas.txt contains the gene names, predicted PASs, and their corresponding fractions in the two conditions.
#' gene.txt contains the genes' PAS number, p values, RED, RLDu, and RLDi scores.
#' @export
#' @import parallel
#' @import GenomicRanges
#' @import GenomicAlignments
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom stats density pchisq fisher.test
#' @importFrom utils write.table
#' @importFrom Rsamtools ScanBamParam
#' @importFrom IRanges subsetByOverlaps IRanges
#' @author Wei Vivian Li, \email{vivian.li@rutgers.edu}
maaper = function(gtf, pas_annotation, output_dir, bam_c1, bam_c2,
                  read_len, ncores = 1,
                  num_pas_thre = 25, frac_pas_thre = 0.05,
                  dist_thre = 600, num_thre = 50,
                  run = "all", subset = NULL, region = "all",
                  gtf_rds = NULL, verbose = FALSE){
  dir.create(output_dir, recursive = T)


  message("Prepare reference genome ...")
  if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/"){
    output_dir = paste0(output_dir, "/")
  }
  if(verbose){
    dir.create((paste0(output_dir, "temp/")), recursive = T)
  }

  if(is.null(gtf_rds)){
    if(run == "all"){
      save_path = paste0(output_dir, "gtf.rds")
      gtf_to_exons_gr(gtf, format = "gtf", save_path)
      gc()
    }
    exons_gr = readRDS(save_path)
  }else{
    exons_gr = readRDS(gtf_rds)
  }

  message("Prepare PAS annotation ...")
  pas_by_gene_single = get_pas_by_gene_single(pas_annotation)
  density_train_path = paste0(output_dir, "density_train.rds")

  save_path = paste0(output_dir, "result.rds")
  wrap(pas_by_gene_single,
       pas_by_gene = pas_annotation, exons_gr,
       bam_c1, bam_c2, density_train_path,
       dist_thre, # distance larger than threshold is not considered in training
       num_thre, # threshold on read number for including in training
       read_len,
       num_pas_thre,
       frac_pas_thre,
       ncores,
       save_path,
       run,
       subset,
       region,
       verbose,
       output_dir)
  gc()
  mapper_write(save_path, output_dir)
}

maaper_predict = function(gtf, pas_annotation, output_dir, bam,
                  read_len, ncores = 1,
                  num_pas_thre = 25, frac_pas_thre = 0.05,
                  dist_thre = 600, num_thre = 50,
                  run = "all", subset = NULL, region = "all",
                  gtf_rds = NULL){
  dir.create(output_dir, recursive = T)

  message("Prepare reference genome ...")
  if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/"){
    output_dir = paste0(output_dir, "/")
  }

  if(is.null(gtf_rds)){
    if(run == "all"){
      save_path = paste0(output_dir, "gtf.rds")
      gtf_to_exons_gr(gtf, format = "gtf", save_path)
      gc()
    }
    exons_gr = readRDS(save_path)
  }else{
    exons_gr = readRDS(gtf_rds)
  }

  message("Prepare PAS annotation ...")
  pas_by_gene_single = get_pas_by_gene_single(pas_annotation)
  density_train_path = paste0(output_dir, "density_train.rds")

  save_path = paste0(output_dir, "result.rds")
  wrap_predict(pas_by_gene_single,
       pas_by_gene = pas_annotation, exons_gr,
       bam,  density_train_path,
       dist_thre, # distance larger than threshold is not considered in training
       num_thre, # threshold on read number for including in training
       read_len,
       num_pas_thre,
       frac_pas_thre,
       ncores,
       save_path,
       run,
       subset,
       region)
  gc()
}
