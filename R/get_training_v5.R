get_pdist_singlePAS = function(bam_path, exons_gr_nover, pas_by_gene_single, 
                               num_thre = 10, dist_thre= 1000, ncores){
  
  genes_single = names(pas_by_gene_single)
  dist_train = mclapply(1:length(pas_by_gene_single), function(i){
    if(i %% 500 == 0) {print(i); gc()}
    gene = genes_single[i]
    type = pas_by_gene_single[[gene]]$Intron.exon.location
    if(type %in% c("Intron"))return(NULL)
    pos = pas_by_gene_single[[gene]]$Pos
    exongr = exons_gr_nover[[gene]]
    str = as.character(strand(exongr)[1])
    chr = as.character(seqnames(exongr)[1])
    ne = length(exongr)
    starts = start(exongr)
    ends = end(exongr)
    exonLens = ends - starts + 1
    if(sum(exonLens) < 400)return(NULL)
    exonLensCum = c(0,cumsum(exonLens))
    
    if(str == "+"){
      gend = max(ends[ne], pos)
      genegr = GRanges(seqnames = chr, 
                       strand = str,
                       ranges = IRanges(start = starts[1], end = gend))
    }
    if(str == "-"){
      gstart = min(starts[1], pos)
      genegr = GRanges(seqnames = chr, 
                       strand = str,
                       ranges = IRanges(start = gstart, end = ends[ne]))
    }
    
    readstarts = get_reads_stranded(gene_gr = genegr, num_thre, bam_path)
    if(is.null(readstarts)) return(NULL)
    
    readgr = GRanges(seqnames = chr, 
                     strand = str,
                     ranges = IRanges(start = readstarts, width = 1))
    readgr = subsetByOverlaps(readgr, exongr)
    readstarts = start(readgr)
    if(str == "+"){
      readstarts =  readstarts[readstarts <= pos]
    }else{
      readstarts =  readstarts[readstarts >= pos]
    }
    if(length(readstarts) <= num_thre) return(NULL)
    if(ne == 1){ dist = abs(pos - readstarts) }
    if(ne > 1){
      #intronLens = starts[-1] - ends[-ne] - 1
      epos = findInterval(pos, starts)
      epos2 = findInterval(pos, ends)
      # maps to exon
      if(epos == epos2 + 1){
        coordOnTxPos = exonLensCum[epos] + pos - starts[epos] + 1
      }else if(epos == ne & epos2 == ne){
        coordOnTxPos = exonLensCum[ne+1] + pos - ends[ne] + 1
      }else if(epos == 0 & epos2 == 0){
        coordOnTxPos = pos - starts[1] + 1
      }else{# not considering intronic APA for now
        return(NULL)
      }
      eread = findInterval(readstarts, starts)
      eread2 = findInterval(readstarts, ends)
      nread = length(readstarts)
      coordOnTx = rep(0, nread)
      # maps to exon
      idx1 = which(eread == eread2 + 1) 
      if(length(idx1) > 0){
        coordOnTx[idx1] = exonLensCum[eread[idx1]] + readstarts[idx1] - starts[eread[idx1]] + 1
      }
      # 3' UTR extension
      idx2 = which(eread == ne & eread2 == ne)
      if(length(idx2) > 0){
        coordOnTx[idx2] = exonLensCum[ne+1] + readstarts[idx2] - ends[ne] + 1
      }
      idx3 = which(eread == 0 & eread2 == 0)
      if(length(idx3) > 0){
        coordOnTx[idx3] = readstarts[idx3] - starts[1] + 1
      }
      # intron (not counting intron length for now)
      idx4 = which(eread == eread2 & eread < ne & eread > 0)
      if(length(idx4) > 0){
        # coordOnTx[idx4] = exonLensCum[eread[idx4]+1] 
        coordOnTx = coordOnTx[-idx4]
      }
      
      dist =  coordOnTxPos - coordOnTx
      if(str == "-"){dist = -dist}
    }
    return(dist)
  }, mc.cores = ncores)
  a = sum(!sapply(dist_train, is.null))
  message(paste(a, "genes used for training ..."))
  
  dist_all = unlist(dist_train)
  b = round(mean(dist_all > dist_thre), digits = 4)
  message(paste(b, "fragments longer than", dist_thre, "..."))
  
  dist_all = dist_all[dist_all <= dist_thre]
  pdist = density(dist_all, from = 0, to = dist_thre, n = dist_thre+1)
  return(pdist)
}


get_reads_prob = function(bam_c1, bam_c2, density_train, conds, exongr, pas_ann, num_thre, dist_thre, read_len, region = "exon-only"){
  
  ## exon stats
  str = pas_ann$Strand[1]
  chr = as.character(seqnames(exongr)[1])
  starts = start(exongr)
  ends = end(exongr)
  ne = length(exongr)
  intronLens = starts[-1] - ends[-ne] - 1
  exonLens = ends - starts + 1
  exonLensCum = c(0,cumsum(exonLens))
  ## pas stats
  pos = pas_ann$Pos
  epos = findInterval(pos, starts)
  epos2 = findInterval(pos, ends)
  npos = length(pos)
  coordOnTxPos = rep(0, npos)
  # maps to exon
  idx1 = which(epos == epos2 + 1) 
  if(length(idx1) > 0){
    coordOnTxPos[idx1] = exonLensCum[epos[idx1]] + pos[idx1] - starts[epos[idx1]] + 1
  }
  # 3' UTR extension
  idx2 = which(epos == ne & epos2 == ne)
  
  if(length(idx2) > 0){
    coordOnTxPos[idx2] = exonLensCum[ne+1] + pos[idx2] - ends[ne] + 1
  }
  idx3 = which(epos == 0 & epos2 == 0)
  if(length(idx3) > 0){
    coordOnTxPos[idx3] = pos[idx3] - starts[1] + 1
  }
  
  # intron 
  idx4 = which(epos == epos2 & epos < ne & epos > 0)
  if(length(idx4) > 0){
    coordOnTxPos[idx4] = pos[idx4] + exonLensCum[ne+1]
  }
  names(coordOnTxPos) = pos
  
  if(str == "+"){
    gend = max(ends[ne], max(pas_ann$Position))
    genegr = GRanges(seqnames = chr, 
                     strand = str,
                     ranges = IRanges(start = starts[1], end = gend))
  }
  if(str == "-"){
    gstart = min(starts[1], min(pas_ann$Position))
    genegr = GRanges(seqnames = chr, 
                     strand = str,
                     ranges = IRanges(start = gstart, end = ends[ne]))
  }
  
  dist_by_pas = lapply(conds, function(con){
    if(con == conds[1]){bam_paths = bam_c1}
    if(con == conds[2]){bam_paths = bam_c2}
    ss = length(bam_paths)
    distmat = lapply(1:ss, function(k){
      bam_path = bam_paths[k]
      readstarts = get_reads_stranded(genegr, num_thre, bam_path)
      if(length(readstarts) <= num_thre) return(NULL)
      
      if(region == "exon-only"){
        readgr = GRanges(seqnames = chr, 
                         strand = str,
                         ranges = IRanges(start = readstarts, width = 1))
        readgr = subsetByOverlaps(readgr, exongr)
        readstarts = start(readgr)
        if(length(readstarts) <= num_thre) return(NULL)
      }
      #if(str == "+"){readstarts = readstarts[readstarts>=starts[1]]}
      #if(str == "-"){readstarts = readstarts[readstarts<=ends[ne]]}
      if(length(readstarts) <= num_thre) return(NULL)
      nread = length(readstarts)
      eread = findInterval(readstarts, starts)
      eread2 = findInterval(readstarts, ends)
      coordOnTx = rep(0, nread)
      # maps to exon
      idx1 = which(eread == eread2 + 1) 
      if(length(idx1) > 0){
        coordOnTx[idx1] = exonLensCum[eread[idx1]] + readstarts[idx1] - starts[eread[idx1]] + 1
      }
      # 3' UTR extension
      idx2 = which(eread == ne & eread2 == ne)
      if(length(idx2) > 0){
        coordOnTx[idx2] = exonLensCum[ne+1] + readstarts[idx2] - ends[ne] + 1
      }
      idx3 = which(eread == 0 & eread2 == 0)
      if(length(idx3) > 0){
        coordOnTx[idx3] = readstarts[idx3] - starts[1] + 1
      }
      # intron 
      idx4 = which(eread == eread2 & eread < ne & eread > 0)
      if(length(idx4) > 0){
        coordOnTx[idx4] = readstarts[idx4] + exonLensCum[ne+1]
      }
      # coordOnTx = coordOnTx[coordOnTx > 0]

      dist = matrix(rep(coordOnTxPos,nread), ncol=npos, byrow=T) -
        matrix(rep(coordOnTx,npos), ncol=npos, byrow=F) 
      if(str == "-"){dist = -dist}
      colnames(dist) = pas_ann$PAS_ID
      return(dist)
    })
    return(distmat)
  })
  names(dist_by_pas) = conds
  ### for filtering
  dist_by_pas_merge = lapply(conds, function(con){
    mat = Reduce(rbind, dist_by_pas[[con]])
    return(mat)
  })
  if(is.null(dist_by_pas_merge[[1]]) | is.null(dist_by_pas_merge[[2]])) return(NULL)
  
  z = min(200, read_len*2)
  tp1 = colSums(dist_by_pas_merge[[1]] <= z & dist_by_pas_merge[[1]] >= 0)
  pas1 = names(tp1)[tp1 > 0]
  tp2 = colSums(dist_by_pas_merge[[2]] <= z & dist_by_pas_merge[[2]] >= 0)
  pas2 = names(tp2)[tp2 > 0]
  pasm = union(pas1, pas2)
  #if(length(pasm) < 2)return(NULL)

  prob_reads_by_pas = lapply(conds, function(con){
    # print(con)
    if(con == conds[1]){bam_paths = bam_c1}
    if(con == conds[2]){bam_paths = bam_c2}
    ss = length(bam_paths)
    probmat = lapply(1:ss, function(k){
      # print(k)
      if(is.null(dist_by_pas[[con]][[k]]))return(NULL)
      distp = density_train[[con]][[k]]$y
      dist = dist_by_pas[[con]][[k]][, pasm, drop = FALSE]
      dval = as.vector(dist)
      pval = rep(0, length(dval))
      pval[dval >= 0 & dval <= dist_thre] = distp[1+dval[dval >= 0 & dval <= dist_thre]]
      pval[dval > dist_thre] = min(distp)/exp(0.05*(dval[dval > dist_thre]-dist_thre))
      pmat = matrix(pval, ncol = ncol(dist))
      colnames(pmat) = colnames(dist)
      pmat = pmat[rowSums(pmat) > 1e-5, , drop = FALSE]
      return(pmat)
    })
    return(probmat)
  })
  names(prob_reads_by_pas) = conds
  nreads = lapply(conds, function(con){
    sapply(prob_reads_by_pas[[con]], function(x){
      if(is.null(x))return(0)
      return(nrow(x))
    })
  })
  if(min(sum(nreads[[1]]), sum(nreads[[2]])) < num_thre)return(NULL)
  names(nreads) = conds
  prob_reads_by_pas = lapply(conds, function(con){
    mat = Reduce(rbind, prob_reads_by_pas[[con]])
    return(mat)
  })
  names(prob_reads_by_pas) = conds
  pas_filtered = lapply(conds, function(con){
    mat = prob_reads_by_pas[[con]]
    colnames(mat)[colSums(mat) > 1e-3]
  })
  pas_filtered = Reduce(union, pas_filtered)
  # if(length(pas_filtered) < 2) return(NULL)
  pas_filtered = sort(pas_filtered)
  prob_reads_by_pas = lapply(prob_reads_by_pas, function(mat){
    mat = mat[, pas_filtered, drop = FALSE]
    return(mat)
  })
  gc()
  return(list(prob = prob_reads_by_pas, nreads = nreads,
              coordOnTxPos = coordOnTxPos))
}


get_reads_prob_predict = function(bam, density_train, exongr, pas_ann, num_thre, dist_thre, read_len, region = "exon-only"){
  ## exon stats
  str = pas_ann$Strand[1]
  chr = as.character(seqnames(exongr)[1])
  starts = start(exongr)
  ends = end(exongr)
  ne = length(exongr)
  intronLens = starts[-1] - ends[-ne] - 1
  exonLens = ends - starts + 1
  exonLensCum = c(0,cumsum(exonLens))
  ## pas stats
  pos = pas_ann$Pos
  epos = findInterval(pos, starts)
  epos2 = findInterval(pos, ends)
  npos = length(pos)
  coordOnTxPos = rep(0, npos)
  # maps to exon
  idx1 = which(epos == epos2 + 1) 
  if(length(idx1) > 0){
    coordOnTxPos[idx1] = exonLensCum[epos[idx1]] + pos[idx1] - starts[epos[idx1]] + 1
  }
  # 3' UTR extension
  idx2 = which(epos == ne & epos2 == ne)
  
  if(length(idx2) > 0){
    coordOnTxPos[idx2] = exonLensCum[ne+1] + pos[idx2] - ends[ne] + 1
  }
  idx3 = which(epos == 0 & epos2 == 0)
  if(length(idx3) > 0){
    coordOnTxPos[idx3] = pos[idx3] - starts[1] + 1
  }
  
  # intron 
  idx4 = which(epos == epos2 & epos < ne & epos > 0)
  if(length(idx4) > 0){
    coordOnTxPos[idx4] = pos[idx4] + exonLensCum[ne+1]
  }
  names(coordOnTxPos) = pos
  
  if(str == "+"){
    gend = max(ends[ne], max(pas_ann$Position))
    genegr = GRanges(seqnames = chr, 
                     strand = str,
                     ranges = IRanges(start = starts[1], end = gend))
  }
  if(str == "-"){
    gstart = min(starts[1], min(pas_ann$Position))
    genegr = GRanges(seqnames = chr, 
                     strand = str,
                     ranges = IRanges(start = gstart, end = ends[ne]))
  }
  
  ss = length(bam)
  dist_by_pas = lapply(1:ss, function(k){
    bam_path = bam[k]
    readstarts = get_reads_stranded(genegr, num_thre, bam_path)
    if(length(readstarts) <= num_thre) return(NULL)
    
    if(region == "exon-only"){
      readgr = GRanges(seqnames = chr, 
                       strand = str,
                       ranges = IRanges(start = readstarts, width = 1))
      readgr = subsetByOverlaps(readgr, exongr)
      readstarts = start(readgr)
      if(length(readstarts) <= num_thre) return(NULL)
    }
    #if(str == "+"){readstarts = readstarts[readstarts>=starts[1]]}
    #if(str == "-"){readstarts = readstarts[readstarts<=ends[ne]]}
    if(length(readstarts) <= num_thre) return(NULL)
    nread = length(readstarts)
    eread = findInterval(readstarts, starts)
    eread2 = findInterval(readstarts, ends)
    coordOnTx = rep(0, nread)
    # maps to exon
    idx1 = which(eread == eread2 + 1) 
    if(length(idx1) > 0){
      coordOnTx[idx1] = exonLensCum[eread[idx1]] + readstarts[idx1] - starts[eread[idx1]] + 1
    }
    # 3' UTR extension
    idx2 = which(eread == ne & eread2 == ne)
    if(length(idx2) > 0){
      coordOnTx[idx2] = exonLensCum[ne+1] + readstarts[idx2] - ends[ne] + 1
    }
    idx3 = which(eread == 0 & eread2 == 0)
    if(length(idx3) > 0){
      coordOnTx[idx3] = readstarts[idx3] - starts[1] + 1
    }
    # intron 
    idx4 = which(eread == eread2 & eread < ne & eread > 0)
    if(length(idx4) > 0){
      coordOnTx[idx4] = readstarts[idx4] + exonLensCum[ne+1]
    }
    # coordOnTx = coordOnTx[coordOnTx > 0]
    
    dist = matrix(rep(coordOnTxPos,nread), ncol=npos, byrow=T) -
      matrix(rep(coordOnTx,npos), ncol=npos, byrow=F) 
    if(str == "-"){dist = -dist}
    colnames(dist) = pas_ann$PAS_ID
    return(dist)
  })

  ### for filtering
  dist_by_pas_merge = Reduce(rbind, dist_by_pas)
  if(is.null(dist_by_pas_merge)) return(NULL)
  
  z = min(200, read_len*2)
  tp1 = colSums(dist_by_pas_merge <= z & dist_by_pas_merge >= 0)
  pasm = names(tp1)[tp1 > 0]
  
  prob_reads_by_pas = lapply(1:ss, function(k){
    # print(k)
    if(is.null(dist_by_pas[[k]]))return(NULL)
    distp = density_train[[k]]$y
    dist = dist_by_pas[[k]][, pasm, drop = FALSE]
    dval = as.vector(dist)
    pval = rep(0, length(dval))
    pval[dval >= 0 & dval <= dist_thre] = distp[1+dval[dval >= 0 & dval <= dist_thre]]
    pval[dval > dist_thre] = min(distp)/exp(0.05*(dval[dval > dist_thre]-dist_thre))
    pmat = matrix(pval, ncol = ncol(dist))
    colnames(pmat) = colnames(dist)
    pmat = pmat[rowSums(pmat) > 1e-5, , drop = FALSE]
    return(pmat)
  })

  nreads = sapply(prob_reads_by_pas, function(x){
    if(is.null(x))return(0)
    return(nrow(x))
  })
  if(sum(nreads) < num_thre) return(NULL)
  
  prob_reads_by_pas = Reduce(rbind, prob_reads_by_pas)

  pas_filtered = colnames(prob_reads_by_pas)[colSums(prob_reads_by_pas) > 1e-3]
  pas_filtered = sort(pas_filtered)
  prob_reads_by_pas = prob_reads_by_pas[ ,pas_filtered, drop = FALSE]
  gc()
  return(list(prob = prob_reads_by_pas, nreads = nreads,
              coordOnTxPos = coordOnTxPos))
}
