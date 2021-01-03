### main function -------------------------
main = function(gene, exongr, pas_ann, bam_c1, bam_c2, density_train, read_len, mode,
                num_thre, dist_thre, sig_fward = 1e-10,
                num_pas_thre, frac_pas_thre, region = "all"){
  conds = c("c1", "c2")
  nExon = length(exongr)
  # if(nrow(pas_ann) == 1) return(NULL)
  
  reads_prob = get_reads_prob(bam_c1, bam_c2, density_train, conds, exongr, pas_ann, num_thre, dist_thre, read_len, region)
  if(is.null(reads_prob)) return(NULL)
  #print("pass1")
  nreads = reads_prob$nreads
  plist = reads_prob$prob
  pasids = colnames(plist[[1]])
  paspos = as.numeric(sapply(strsplit(pasids, split = "\\:"), function(x) x[2]))
  npas = length(colnames(plist[[1]]))
  
  if(npas == 1){
    pas.type = pas_ann[pas_ann$PAS_ID == pasids, "Intron.exon.location"]
    res = list(gene = gene, pas = pasids, pas.type = pas.type,
               nreads = nreads, npas = 1)
    return(res)
  }
  
  alp = matrix(1, nrow = 1, ncol = 2)
  loglik_init = sapply(1:npas, function(j){
    return( get_loglik(plist[[1]][,j,drop = FALSE], plist[[2]][,j,drop = FALSE],  alp_hat = alp) )
  })
  loglik_null = max(loglik_init)
  idx_pas = which.max(loglik_init)
  select_pas = select_forward(max_step = 2*nExon, idx_init = idx_pas, 
                              pas_pos = paspos, plist,
                              loglik = loglik_null, sig_fward, mode)
  if(is.null(select_pas))return(NULL)
  #print("pass2")
  idx_pas = select_pas$idx_sel
  alpha = select_pas$alp
  
  ind1 = (max(nreads$c1) * alpha[,1] >= num_pas_thre) |
    (max(nreads$c2) * alpha[,2] >= num_pas_thre)
  ind2 = (alpha[,1] >= frac_pas_thre) | (alpha[,2] >= frac_pas_thre)
  ind = which(ind1 & ind2)
  if(length(ind) == 0) return(NULL)
  idx_pas = idx_pas[ind]
  if(length(ind) == 1){
    pas.type = pas_ann$Intron.exon.location[match(pasids[idx_pas], pas_ann$PAS_ID)]
    res = list(gene = gene, pas = pasids[idx_pas], pas.type = pas.type,
               nreads = nreads, npas = 1)
    return(res)
  }
  ### alpha and loglik assuming DE
  de_alt = get_de_alt(idx_sel = idx_pas, plist)
  loglik_de_alt = de_alt$loglik
  alp_alt = de_alt$alp
  rownames(alp_alt) = pasids[idx_pas]
  de_null = get_de_null(idx_sel = idx_pas, plist, iter = 500)
  alp_null = de_null$alp
  loglik_de_null = de_null$loglik
  Q = 2 * (loglik_de_alt - loglik_de_null)
  pval = pchisq(Q, df = length(idx_pas), lower.tail = F)
  # TV = 0.5*sum(abs(alp_alt[,1] - alp_alt[,2]))
  pasids = names(alp_null)
  pas.type = pas_ann$Intron.exon.location[match(pasids, pas_ann$PAS_ID)]
  
  res = list(gene = gene, alp_alt = alp_alt, 
             pval = pval, nreads = nreads, 
             pas.type = pas.type, npas = length(pas.type))

  ### supplemental info
  frac = res$alp_alt
  if(is.null(frac)){return(res)}
  
  pas.type = pas.type[order(rownames(frac))]
  frac = frac[order(rownames(frac)), ]
  
  if(all(pas.type == "Single exon")){
    RLD = get_RLD(exongr, frac)
  }else if(sum(pas.type == "3' most exon") >= 2){
    frac.rld = frac[pas.type == "3' most exon", ]
    RLD = get_RLD(exongr, frac.rld)
  }else{RLD = NA}
  res$RLD.utr = RLD
  
  if(sum(pas.type != "3' most exon") >= 1 & sum(pas.type == "3' most exon") >= 1){
    RLD = get_RID(exongr, frac, pas.type)
  }else{RLD = NA}
  res$RLD.all = RLD
  
  if(nrow(frac) == 2){
    #top2.exp = rownames(frac)
    top2.diff = rownames(frac)
  }else{
    # top2.exp = rownames(frac)[order(rowSums(frac), decreasing = T)[1:2]]
    # top2.exp = sort(top2.exp)
    top2.diff = rownames(frac)[order(abs(frac[,1]-frac[,2]), decreasing = T)[1:2]]
    top2.diff = sort(top2.diff)
  }
  ### calculate RED using most 2 abundant PASs
  # frac2 = frac[top2.exp, ]
  # nreads = sapply(res$nreads, sum)
  # nr = ceiling(sweep(frac2, 2, nreads, "*"))
  # ftest = fisher.test(t(nr))
  # pval.fisher = ftest$p.value
  # RED = log2(frac2[2,2]/frac2[1,2]) - log2(frac2[2,1]/frac2[1,1])
  # if(pas_ann$Strand[1] == "-"){RED = -RED}
  # names(RED) = NULL
  # exp2 = list(RED = RED, pval = pval.fisher, top2 = top2.exp)
  # res$exp2 = exp2
  ### calculate RED using most 2 differential PASs
  frac2 = frac[top2.diff, ]
  nreads = sapply(res$nreads, sum)
  nr = ceiling(sweep(frac2, 2, nreads, "*"))
  ftest = fisher.test(t(nr))
  pval.fisher = ftest$p.value
  RED = log2(frac2[2,2]/frac2[1,2]) - log2(frac2[2,1]/frac2[1,1])
  if(pas_ann$Strand[1] == "-"){RED = -RED}
  names(RED) = NULL
  diff2 = list(RED = RED, pval = pval.fisher, top2 = top2.diff)
  res$diff2 = diff2
  # if(all(top2.exp == top2.diff)){
  #   res$diff2 = exp2
  # }else{
  # }
  return(res)
}


### wrap function -------------------------
wrap = function(pas_by_gene_single, pas_by_gene, exons_gr,
                bam_c1, bam_c2, density_train_path,
                dist_thre = 600, # distance larger than threshold is not considered in training
                num_thre = 10, # threshold on read number for including in training
                read_len = 76,
                num_pas_thre = 10,
                frac_pas_thre = 0.01,
                ncores = 1,
                save_path,
                run = "all",
                subset = NULL,
                region = "all",
                verbose = FALSE,
                output_dir = NULL
){
  ### PAS annotation ---------------------
  GENEID = sapply(pas_by_gene, function(x) x$Ensemble.ID[1])
  GENESB = names(GENEID)
  ### reference genome -------------------
  geneid = names(exons_gr)
  geneid = sapply(strsplit(geneid, split = "\\."), function(x) x[1])
  genesb = GENESB[match(geneid, GENEID)]
  exons_gr = exons_gr[!is.na(genesb)]
  genesb = genesb[!is.na(genesb)]
  message(paste(length(genesb), "genes!"))
  names(exons_gr) = genesb
  
  pas_by_gene_single = pas_by_gene_single[names(pas_by_gene_single) %in% genesb]
  pas_by_gene = pas_by_gene[names(pas_by_gene) %in% genesb]
  
  ### training --------------------------
  conds = c("c1", "c2")
  if(run == "all"){
    density_train = lapply(conds, function(con){
      if(con == "c1"){bam_paths = bam_c1}
      if(con == "c2"){bam_paths = bam_c2}
      ss = length(bam_paths)
      density_con = lapply(1:ss, function(k){
        message(paste("Start training for condition", con, "-", "sample", k,  "..."))
        bam_path = bam_paths[k]
        pdist = get_pdist_singlePAS(bam_path, exons_gr, pas_by_gene_single, num_thre, dist_thre, ncores)
        gc()
        return(pdist)
      })
    })
    names(density_train) = conds
    saveRDS(density_train, file = density_train_path)
    gc()
  }
  if(run == "skip-train"){
    density_train = readRDS(density_train_path)
  }
  mode = lapply(conds, function(con){
    if(con == "c1"){bam_paths = bam_c1}
    if(con == "c2"){bam_paths = bam_c2}
    ss = length(bam_paths)
    mod = sapply(1:ss, function(k){
      dens = density_train[[con]][[k]]
      dens$x[which.max(dens$y)]
    })
    return(mod)
  })
  
  mode = ceiling(mean(unlist(mode)))
  
  message("Calculate APA events ...")
  if(is.null(subset)){usegenes = GENESB
  }else{
    message(paste(length(subset), "genes in subset!"))
    usegenes = intersect(subset,GENESB)
    message(paste(length(usegenes), "genes in subset and data!"))
    if(length(usegenes) == 0) message("Stop!")
  }
  
  result = mclapply(usegenes, function(gene){
    print(gene)
    exongr = exons_gr[[gene]]
    if(is.null(exongr)) return(NULL)
    pas_ann = pas_by_gene[[gene]]
    str = pas_ann$Strand[1]
    if(str == "+"){
      idx = which(pas_ann$Position > start(exongr)[1])
      if(length(idx) == 0){return(NULL)}
      pas_ann = pas_ann[idx, , drop = FALSE]
    }
    
    if(str == "-"){
      idx = which(pas_ann$Position < max(end(exongr)))
      if(length(idx) == 0){return(NULL)}
      pas_ann = pas_ann[idx, , drop = FALSE]
    }
    
    res = try(main(gene, exongr, pas_ann, bam_c1, bam_c2, density_train, read_len, mode,
                   num_thre, dist_thre, sig_fward = 1e-10,
                   num_pas_thre, frac_pas_thre, region), silent = TRUE)
    gc()
    if(verbose) saveRDS(res, file = paste0(output_dir, "temp/", gene, ".rds"))
    return(res)
  }, mc.cores = ncores)
  names(result) = usegenes
  saveRDS(result, file = save_path)
  return(0)
}

### write function -------------------------
mapper_write = function(save_path, output_dir){
  result = readRDS(save_path)
  ind = which(sapply(result, function(x) class(x) != "try-error"))
  result = result[ind]
  result = result[!sapply(result, is.null)]
  message(paste("Writing results for", length(result), "genes ..."))
          
  npas = sapply(result, function(x) x$npas)
  ### multiple pas
  res1 = result[npas > 1]
  da.gene = lapply(1:length(res1), function(i){
    res = res1[[i]]
    da = data.frame(gene = res$gene, npas = res$npas, 
                    pval = signif(res$pval, digits = 3), 
                    RLD.u = signif(res$RLD.utr, digits = 3), 
                    RLD.i = signif(res$RLD.all, digits = 3), 
                    RED = signif(res$diff2$RED, digits = 3), 
                    pval.fisher = signif(res$diff2$pval, digits = 3),
                    RED.gene1 = res$diff2$top2[1], RED.gene2 = res$diff2$top2[2],
                    t(res$nreads$c1), t(res$nreads$c2))
    return(da)
  })
  da.gene = Reduce(rbind, da.gene)
  l = (ncol(da.gene)-9)/2
  colnames(da.gene) = c(colnames(da.gene)[1:9],
                       paste0("nread.c1.", 1:l),
                       paste0("nread.c2.", 1:l))
  write.table(da.gene, paste0(output_dir, "gene.txt"),
              sep="\t", row.names=FALSE, quote = F)
  
  da.pas = lapply(1:length(res1), function(i){
    res = res1[[i]]
    count = data.frame(gene = res$gene, pas = rownames(res$alp_alt), type = res$pas.type, 
                       frac.c1 = signif(res$alp_alt[,1], digits = 3), 
                       frac.c2 = signif(res$alp_alt[,2], digits= 3))
    return(count)
  })
  da.pas = Reduce(rbind, da.pas)
  write.table(da.pas, paste0(output_dir, "pas.txt"),
              sep="\t", row.names=FALSE, quote = F)
  
  ### single pas
  res1 = result[npas == 1]
  da.gene = lapply(1:length(res1), function(i){
    res = res1[[i]]
    da = data.frame(gene = res$gene, npas = 1, pval = NA, 
                    RLD.u = NA, RLD.i = NA, RED = NA, pval.fisher = NA,
                    RED.gene1 = NA, RED.gene2 = NA,
                    t(res$nreads$c1), t(res$nreads$c2))
    return(da)
  })
  da.gene = Reduce(rbind, da.gene)
  write.table(da.gene, paste0(output_dir, "gene.txt"),
              sep="\t", row.names=FALSE, col.names = FALSE,
              quote = F, append = T)
  
  da.pas = lapply(1:length(res1), function(i){
    res = res1[[i]]
    count = data.frame(res$gene, res$pas, res$pas.type, 1, 1)
    return(count)
  })
  da.pas = Reduce(rbind, da.pas)
  # colnames(da.pas) = c("gene", "pas", "type", "frac.c1", "frac.c2",
  #                      paste0("nread.c1.", 1:length(res$nreads$c1)),
  #                      paste0("nread.c2.", 1:length(res$nreads$c2)))
  write.table(da.pas, paste0(output_dir, "pas.txt"),
              sep="\t", row.names=FALSE, col.names = FALSE,
              quote = F, append = T)
  return("Done!")
}

