
get_reads_stranded = function(gene_gr, num_thre = 10, bam_path){
  param = ScanBamParam(which = gene_gr)
  str = as.character(strand(gene_gr))[1]
  ### Given single-end data
  reads = readGAlignments(bam_path, param = param)
  if(length(reads) < num_thre) return(NULL)
  reads = reads[as.character(strand(reads)) == str]
  if(length(reads) < num_thre) return(NULL)

  if(str == "+"){
    reads_start = start(reads)
  }
  if(str == "-"){
    reads_start = end(reads)
  }

  return(reads_start)
}

