get_pas_by_gene_single = function(pas_by_gene){
  pas_by_gene_single = pas_by_gene[sapply(pas_by_gene, nrow) == 1]
  pas_by_gene_single = pas_by_gene_single[sapply(pas_by_gene_single, function(x) x$LOCATION) != "Intron"]
  length(pas_by_gene_single)
  names(pas_by_gene_single) = sapply(pas_by_gene_single, function(x) x$Gene.Symbol[1])
  return(pas_by_gene_single)
}