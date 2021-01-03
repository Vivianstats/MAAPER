

gtf_to_exons_gr = function(gtf_path, format = format, save_path){
  txdb = makeTxDbFromGFF(gtf_path, format)
  txdb = keepStandardChromosomes(txdb)
  # dropSeqlevels(txdb, "chrM")
  exons_list_by_gene = exonsBy(txdb, by="gene")
  exons_list_by_gene = reduce(exons_list_by_gene)
  saveRDS(exons_list_by_gene, file = save_path)
  return(0)
}

