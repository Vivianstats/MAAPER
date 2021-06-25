get_pas_by_gene_single = function(pas_by_gene){
  pas_by_gene_single = pas_by_gene[sapply(pas_by_gene, nrow) == 1]
  pas_by_gene_single = pas_by_gene_single[sapply(pas_by_gene_single, function(x) x$LOCATION) != "Intron"]
  length(pas_by_gene_single)
  names(pas_by_gene_single) = sapply(pas_by_gene_single, function(x) x$Gene.Symbol[1])
  return(pas_by_gene_single)
}

trunc_matrix = function(Smat){
  eSmat = eigen(Smat)
  eSmat$values[eSmat$values < 0] = 0
  Smat = eSmat$vectors %*% diag(eSmat$values) %*% t(eSmat$vectors)
  return(Smat)
}

paired_test = function(X1, X2){
  npas = nrow(X1)
  n = ncol(X1)

  N1 = colSums(X1)
  N2 = colSums(X2)
  NN1 = sum(N1)
  NN2 = sum(N2)

  mi1 = t(t(X1) / N1)
  mi2 = t(t(X2) / N2)
  fraclist = lapply(1:n, function(i){
    cbind(mi1[,i], mi2[, i])
  })

  if(NN1 <= n | NN2 <= n){
    return(list(pvalue = NA, fraclist = fraclist))
  }else{
    idx = which(N1 > 0 & N2 > 0)
    if(length(idx) < n){
      if(length(idx) == 1) return(list(pvalue = NA, fraclist = fraclist))
      X1 = X1[, idx]
      X2 = X2[, idx]
      n = ncol(X1)
    }
    if(n < 3){return(list(pvalue = NA, fraclist = fraclist))}
    if(n - npas < 1){
      a = rowSums(X1) + rowSums(X2)
      idx = order(a, decreasing = T)[1:(n-1)]
      X1 = X1[idx, ,drop=FALSE]
      X2 = X2[idx, ,drop=FALSE]
      npas = nrow(X1)
    }

    N1 = colSums(X1)
    N2 = colSums(X2)
    NN1 = sum(N1)
    NN2 = sum(N2)

    mi1 = t(t(X1) / N1)
    mi2 = t(t(X2) / N2)

    Mc1 = (NN1 - sum(N1^2) / NN1) / (n-1)
    Mc2 = (NN2 - sum(N2^2) / NN2) / (n-1)

    m1 = rowSums(X1) / NN1
    m2 = rowSums(X2) / NN2
    mi1 = t(t(X1) / N1)
    mi2 = t(t(X2) / N2)
    diff1 = mi1 - m1
    diff2 = mi2 - m2
    S1 = diff1 %*% diag(N1) %*% t(diff1) / (n-1)
    S2 = diff2 %*% diag(N2) %*% t(diff2) / (n-1)
    G1 = (diag(mi1 %*% N1) - mi1 %*% diag(N1) %*% t(mi1)) / (NN1-n)
    G2 = (diag(mi2 %*% N2) - mi2 %*% diag(N2) %*% t(mi2)) / (NN2-n)
    V = diff1 %*% diag(N1+N2) %*% t(diff2) /(Mc1+Mc2) / (n-1)

    V = (V + t(V)) * sum(N1*N2) / (NN1*NN2)
    Sig1 = S1 * sum(N1^2) / Mc1 / (NN1^2) + G1 * (Mc1 - sum(N1^2) / NN1) / NN1 / Mc1
    Sig2 = S2 * sum(N2^2) / Mc2 / (NN2^2) + G2 * (Mc2 - sum(N2^2) / NN2) / NN2 / Mc2
    Sig1 = trunc_matrix(Sig1)
    Sig2 = trunc_matrix(Sig2)

    Smat = Sig1 + Sig2 - V
    Smat = trunc_matrix(Smat)

    stat = t(m1-m2) %*% ginv(Smat) %*% (m1-m2) * (n-npas+1)/(n-1)/(npas-1)
    pval = 1-pf(as.numeric(stat), npas-1, n-npas+1)
  }
  return(list(Smat = Smat, pvalue = pval, fraclist = fraclist))
}

write_bed = function(save_path, output_dir, read_len, paired, ns){
  result = readRDS(save_path)
  ind = which(sapply(result, function(x) class(x) != "try-error"))
  result = result[ind]
  result = result[!sapply(result, is.null)]
  message(paste("Writing BED for", length(result), "genes ..."))

  if(paired == FALSE){
    dabed = lapply(1:length(result), function(i){
      res = result[[i]]
      if(res$npas > 1){
        pasids = rownames(res$alp_alt)
        pasids = strsplit(pasids, split = "\\:")
        chr = sapply(pasids, function(x) x[1])
        str = pasids[[1]][3]
        if(str == "+"){
          chrEnd = as.numeric(sapply(pasids, function(x) x[2]))
          chrStart = chrEnd - read_len + 1
        }else{
          chrStart = as.numeric(sapply(pasids, function(x) x[2]))
          chrEnd = chrStart + read_len - 1
        }
        #count_c1 = round(sum(res$nreads$c1) * res$alp_alt[,1])
        #count_c2 = round(sum(res$nreads$c2) * res$alp_alt[,2])
        count_c1 = res$alp_alt[,1]
        count_c2 = res$alp_alt[,2]
        da_c1 = data.frame(chr, chrStart, chrEnd, count_c1)
        da_c2 = data.frame(chr, chrStart, chrEnd, count_c2)
      }else{
        pasids = res$pas
        pasids = strsplit(pasids, split = "\\:")
        chr = sapply(pasids, function(x) x[1])
        str = pasids[[1]][3]
        if(str == "+"){
          chrEnd = as.numeric(sapply(pasids, function(x) x[2]))
          chrStart = chrEnd - read_len + 1
        }else{
          chrStart = as.numeric(sapply(pasids, function(x) x[2]))
          chrEnd = chrStart + read_len - 1
        }
        # count_c1 = sum(res$nreads$c1)
        # count_c2 = sum(res$nreads$c2)
        count_c1 = 1
        count_c2 = 1
        da_c1 = data.frame(chr, chrStart, chrEnd, count_c1)
        da_c2 = data.frame(chr, chrStart, chrEnd, count_c2)
      }
      return(list(c1 = da_c1, c2 = da_c2, str = str))
    })
    idx1 = sapply(dabed, function(x) x$str == "+")
    idx2 = sapply(dabed, function(x) x$str == "-")

    dabed_c1 = Reduce(rbind, lapply(dabed[idx1], function(x) x$c1))
    write.table(dabed_c1, paste0(output_dir, "c1_forward.bedgraph"),
                sep="\t", row.names=FALSE, col.names=FALSE, quote = F)
    dabed_c2 = Reduce(rbind, lapply(dabed[idx1], function(x) x$c2))
    write.table(dabed_c2, paste0(output_dir, "c2_forward.bedgraph"),
                sep="\t", row.names=FALSE, col.names=FALSE, quote = F)

    dabed_c1 = Reduce(rbind, lapply(dabed[idx2], function(x) x$c1))
    write.table(dabed_c1, paste0(output_dir, "c1_reverse.bedgraph"),
                sep="\t", row.names=FALSE, col.names=FALSE, quote = F)
    dabed_c2 = Reduce(rbind, lapply(dabed[idx2], function(x) x$c2))
    write.table(dabed_c2, paste0(output_dir, "c2_reverse.bedgraph"),
                sep="\t", row.names=FALSE, col.names=FALSE, quote = F)
  }

  if(paired == TRUE){
    dabed = lapply(1:length(result), function(i){
      res = result[[i]]
      if(res$npas > 1){
        pasids = rownames(res$alp_alt)
        pasids = strsplit(pasids, split = "\\:")
        chr = sapply(pasids, function(x) x[1])
        str = pasids[[1]][3]
        if(str == "+"){
          chrEnd = as.numeric(sapply(pasids, function(x) x[2]))
          chrStart = chrEnd - read_len + 1
        }else{
          chrStart = as.numeric(sapply(pasids, function(x) x[2]))
          chrEnd = chrStart + read_len - 1
        }
        da_by_samples = lapply(1:ns, function(k){
          #count_c1 = round(sum(res$nreads$c1) * res$alp_alt[,1])
          #count_c2 = round(sum(res$nreads$c2) * res$alp_alt[,2])
          count_c1 = res$fraclist[[k]][,1]
          count_c2 = res$fraclist[[k]][,2]
          da_c1 = data.frame(chr, chrStart, chrEnd, count_c1)
          da_c2 = data.frame(chr, chrStart, chrEnd, count_c2)
          return(list(c1 = da_c1, c2 = da_c2))
        })
      }else{
        pasids = res$pas
        pasids = strsplit(pasids, split = "\\:")
        chr = sapply(pasids, function(x) x[1])
        str = pasids[[1]][3]
        if(str == "+"){
          chrEnd = as.numeric(sapply(pasids, function(x) x[2]))
          chrStart = chrEnd - read_len + 1
        }else{
          chrStart = as.numeric(sapply(pasids, function(x) x[2]))
          chrEnd = chrStart + read_len - 1
        }
        da_by_samples = lapply(1:ns, function(k){
          # count_c1 = sum(res$nreads$c1)
          # count_c2 = sum(res$nreads$c2)
          count_c1 = 1
          count_c2 = 1
          da_c1 = data.frame(chr, chrStart, chrEnd, count_c1)
          da_c2 = data.frame(chr, chrStart, chrEnd, count_c2)
          return(list(c1 = da_c1, c2 = da_c2))
        })
      }
      return(list(da_by_samples = da_by_samples, str = str))
    })
    gc()
    idx1 = sapply(dabed, function(x) x$str == "+")
    idx2 = sapply(dabed, function(x) x$str == "-")
    for(k in 1:ns){
      dabed_c1 = Reduce(rbind, lapply(dabed[idx1], function(x) x$da_by_samples[[k]]$c1))
      dabed_c1 = dabed_c1[complete.cases(dabed_c1), ]
      write.table(dabed_c1, paste0(output_dir, "c1_sample", k, "_forward.bedgraph"),
                  sep="\t", row.names=FALSE, col.names=FALSE, quote = F)

      dabed_c2 = Reduce(rbind, lapply(dabed[idx1], function(x) x$da_by_samples[[k]]$c2))
      dabed_c2 = dabed_c2[complete.cases(dabed_c2), ]
      write.table(dabed_c2, paste0(output_dir, "c2_sample", k, "_forward.bedgraph"),
                  sep="\t", row.names=FALSE, col.names=FALSE, quote = F)

      dabed_c1 = Reduce(rbind, lapply(dabed[idx2], function(x) x$da_by_samples[[k]]$c1))
      dabed_c1 = dabed_c1[complete.cases(dabed_c1), ]
      write.table(dabed_c1, paste0(output_dir, "c1_sample", k, "_reverse.bedgraph"),
                  sep="\t", row.names=FALSE, col.names=FALSE, quote = F)
      dabed_c2 = Reduce(rbind, lapply(dabed[idx2], function(x) x$da_by_samples[[k]]$c2))
      dabed_c2 = dabed_c2[complete.cases(dabed_c2), ]
      write.table(dabed_c2, paste0(output_dir, "c2_sample", k, "_reverse.bedgraph"),
                  sep="\t", row.names=FALSE, col.names=FALSE, quote = F)
    }
  }
  return(0)
}


