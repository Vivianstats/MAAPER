### For prediction
get_loglik_pred = function(pmat1, alp_hat){
  tp1 = pmat1 %*% alp_hat[, 1, drop = FALSE] #c1
  tp1[tp1 <= 1e-100] = 1e-100
  return(sum(log(tp1)))
}
est_alp_forward_pred = function(pmat1, init, iter = 500){
  J = nrow(init)
  alp_crr = init
  ii = 0
  while(ii < iter){
    # print(ii)
    alp_new = alp_crr
    n1 = nrow(pmat1)
    denom1 = pmat1 %*% alp_crr[,1, drop = FALSE]
    denom1 = matrix(rep(denom1, J), ncol = J, byrow = F)
    num1 = pmat1 * matrix(rep(alp_crr[,1], n1), nrow = n1, byrow = T)
    alp_new[,1] = colMeans(num1/denom1, na.rm = T)

    error = mean((alp_new - alp_crr)^2)
    #print(error)
    if(error < 1e-5){
      alp_crr = alp_new
      break;
    }else{alp_crr = alp_new}
    ii = ii + 1
  }
  return(alp_crr)
}

select_forward_pred = function(max_step, idx_init, pas_pos,
                          plist, loglik, sig_fward, mode){
  steps = 0
  JJ = ncol(plist)
  loglik_null = loglik
  idx_sel = idx_init
  idx_out = c()
  ncrr = 1
  while(steps < max_step){
    loglik_forw = sapply(1:JJ, function(k){
      if(k %in% idx_sel)return(-1e200)
      if(k %in% idx_out)return(-1e200)
      idx_crr = c(idx_sel, k)
      init = matrix(1/(ncrr+1), ncol = 2, nrow = ncrr+1)
      pmat1 = plist[, idx_crr, drop = FALSE]
      alp_curr = est_alp_forward_pred(pmat1, init, iter = 500)
      loglik_curr = get_loglik_pred(pmat1, alp_hat = alp_curr) 
      return(loglik_curr)
    })
    loglik_alt = max(loglik_forw)
    idx_new = which.max(loglik_forw)
    
    loglik_diff = loglik_alt - loglik_null
    dist_to_sel = min(abs(pas_pos[idx_new] - pas_pos[idx_sel]))
    if(dist_to_sel < 100){ # < mode
      idx_out = c(idx_out, idx_new)
      if(length(idx_sel) == JJ-length(idx_out)){break}
      next
    }
    Q = 2*loglik_diff
    pval = pchisq(Q, 1, lower.tail = F)
    if(pval <= sig_fward){
      ncrr = ncrr + 1
      idx_sel = c(idx_sel, idx_new)
      loglik_null =  loglik_alt
    }else{break}
    if(length(idx_sel) == JJ-length(idx_out)){break}
    steps = steps + 1
  }
  if(length(idx_sel) < 2)return(list(idx_sel = idx_sel, 
                                     alp = matrix(1, nrow=1, ncol = 2)))
  
  init = matrix(1/ncrr, ncol = 2, nrow = ncrr)
  pmat1 = plist[, idx_sel, drop = FALSE]
  alp_curr = est_alp_forward_pred(pmat1, init, iter = 500)
  return(list(idx_sel = idx_sel, alp = alp_curr))
  
}

### For APA
get_loglik = function(pmat1, pmat2, alp_hat){
  tp1 = pmat1 %*% alp_hat[, 1, drop = FALSE] #c1
  tp2 = pmat2 %*% alp_hat[, 2, drop = FALSE] #c1
  tp1[tp1 <= 1e-100] = 1e-100
  tp2[tp2 <= 1e-100] = 1e-100
  return(sum(log(tp1)) + sum(log(tp2)))
}

est_alp_forward = function(pmat1, pmat2, init, iter = 500){
  J = nrow(init)
  alp_crr = init
  ii = 0
  while(ii < iter){
    # print(ii)
    alp_new = alp_crr
    n1 = nrow(pmat1)
    denom1 = pmat1 %*% alp_crr[,1, drop = FALSE]
    denom1 = matrix(rep(denom1, J), ncol = J, byrow = F)
    num1 = pmat1 * matrix(rep(alp_crr[,1], n1), nrow = n1, byrow = T)
    alp_new[,1] = colMeans(num1/denom1, na.rm = T)
    # alp_new[,1] = colMeans((num1[denom1>0])/(denom1[denom1>0]))
    
    n2 = nrow(pmat2)
    denom2 = pmat2 %*% alp_crr[,2, drop = FALSE]
    denom2 = matrix(rep(denom2, J), ncol = J, byrow = F)
    num2 = pmat2 * matrix(rep(alp_crr[,2], n2), nrow = n2, byrow = T)
    alp_new[,2] = colMeans(num2/denom2, na.rm = T)
    # alp_new[,2] = colMeans((num2[denom2>0])/(denom2[denom2>0]))
    
    error = mean((alp_new - alp_crr)^2)
    #print(error)
    if(error < 1e-5){
      alp_crr = alp_new
      break;
    }else{alp_crr = alp_new}
    ii = ii + 1
  }
  return(alp_crr)
}


select_forward = function(max_step, idx_init, pas_pos,
                          plist, loglik, sig_fward, mode){
  steps = 0
  JJ = ncol(plist[[1]])
  loglik_null = loglik
  idx_sel = idx_init
  idx_out = c()
  ncrr = 1
  while(steps < max_step){
    # print(steps)
    # print(idx_sel)
    # print(idx_out)
    loglik_forw = sapply(1:JJ, function(k){
      if(k %in% idx_sel)return(-1e200)
      if(k %in% idx_out)return(-1e200)
      idx_crr = c(idx_sel, k)
      init = matrix(1/(ncrr+1), ncol = 2, nrow = ncrr+1)
      pmat1 = plist[[1]][, idx_crr, drop = FALSE]
      pmat2 = plist[[2]][, idx_crr, drop = FALSE]
      alp_curr = est_alp_forward(pmat1, pmat2, init, iter = 500)
      loglik_curr = get_loglik(plist[[1]][,idx_crr,drop = FALSE], 
                               plist[[2]][,idx_crr,drop = FALSE],  
                               alp_hat = alp_curr) 
      return(loglik_curr)
    })
    loglik_alt = max(loglik_forw)
    idx_new = which.max(loglik_forw)
    
    loglik_diff = loglik_alt - loglik_null
    dist_to_sel = min(abs(pas_pos[idx_new] - pas_pos[idx_sel]))
    if(dist_to_sel < 100){ # < mode
      idx_out = c(idx_out, idx_new)
      if(length(idx_sel) == JJ-length(idx_out)){break}
      next
    }
    Q = 2*loglik_diff
    pval = pchisq(Q, 1, lower.tail = F)
    if(pval <= sig_fward){
      ncrr = ncrr + 1
      idx_sel = c(idx_sel, idx_new)
      loglik_null =  loglik_alt
      # print(pval)
      # print(loglik_null)
      
      ## for testing
      # init = matrix(1/ncrr, ncol = 2, nrow = ncrr)
      # pmat1 = plist[[1]][, idx_sel]
      # pmat2 = plist[[2]][, idx_sel]
      # alp_curr = est_alp_forward(pmat1, pmat2, init, iter = 500)
      # print(alp_curr)
    }else{break}
    if(length(idx_sel) == JJ-length(idx_out)){break}
    steps = steps + 1
  }
  if(length(idx_sel) < 2)return(list(idx_sel = idx_sel, 
                                     alp = matrix(1, nrow=1, ncol = 2)))
  ## for testing
  init = matrix(1/ncrr, ncol = 2, nrow = ncrr)
  pmat1 = plist[[1]][, idx_sel, drop = FALSE]
  pmat2 = plist[[2]][, idx_sel, drop = FALSE]
  alp_curr = est_alp_forward(pmat1, pmat2, init, iter = 500)
  return(list(idx_sel = idx_sel, alp = alp_curr))

  }

### alpha and loglik assuming DE
get_de_alt = function(idx_sel, plist){
  ncrr = length(idx_sel)
  init = matrix(1/ncrr, ncol = 2, nrow = ncrr)
  pmat1 = plist[[1]][, idx_sel]
  pmat2 = plist[[2]][, idx_sel]
  alp_curr = est_alp_forward(pmat1, pmat2, init, iter = 500)
  loglik_curr = get_loglik(pmat1, pmat2, alp_hat = alp_curr)
  return(list(alp = alp_curr, loglik = loglik_curr))
}

### alpha and loglik assuming DE
get_de_null = function(idx_sel, plist, iter = 500){
  ncrr = length(idx_sel)
  init = matrix(1/ncrr, ncol = 2, nrow = ncrr)
  pmat1 = plist[[1]][, idx_sel]
  pmat2 = plist[[2]][, idx_sel]
  pmat = rbind(pmat1, pmat2)
  ### get alpha
  alp_crr = rep(1/ncrr, ncrr)
  ii = 0
  while(ii < iter){
    # print(ii)
    alp_new = alp_crr
    n = nrow(pmat)
    denom = pmat %*% matrix(alp_crr, ncol = 1)
    denom = matrix(rep(denom, ncrr), ncol = ncrr, byrow = F)
    num = pmat * matrix(rep(alp_crr, n), nrow = n, byrow = T)
    alp_new = colMeans(num/denom, na.rm = T)
    error = mean((alp_new - alp_crr)^2)
    #print(error)
    if(error < 1e-5){
      alp_crr = alp_new
      break;
    }else{alp_crr = alp_new}
    ii = ii + 1
  }
  ### get loglik
  tp1 = pmat %*% matrix(alp_crr, ncol = 1)
  tp1[tp1 <= 1e-100] = 1e-100
  loglik_crr = sum(log(tp1))
  return(list(alp = alp_crr, loglik = loglik_crr))
}


### calculate RLD (relative length difference)
get_RLD = function(exongr, frac){
  str = as.character(strand(exongr)[1])
  
  pos = strsplit(rownames(frac), split="\\:")
  pos = as.numeric(sapply(pos, function(x) x[2]))

  if(str == "+")coord = pos - min(pos) + 1
  if(str == "-")coord = max(pos) - pos + 1
  RLD = log2(sum(frac[,2] * coord)/sum(frac[,1] * coord))
  return(RLD)
}

get_RID = function(exongr, frac, pas.type){
  pos = strsplit(rownames(frac), split="\\:")
  pos = as.numeric(sapply(pos, function(x) x[2]))
  
  ## exon stats
  str = as.character(strand(exongr)[1])
  
  if(str == "+"){
    coord = pos - min(pos) + 1
  }
  if(str == "-"){
    coord = max(pos) - pos + 1
  }
  idxi = which(pas.type != "3' most exon")
  idxu = which(pas.type == "3' most exon")
  ratio2 = sum(frac[idxu, 2]*coord[idxu])/sum(frac[idxi,2]*coord[idxi])
  ratio1 = sum(frac[idxu, 1]*coord[idxu])/sum(frac[idxi,1]*coord[idxi])
  RLD = log2(ratio2/ratio1)
  return(RLD)
}

# get_RGD = function(exongr, frac){
#   pos = strsplit(rownames(frac), split="\\:")
#   pos = as.numeric(sapply(pos, function(x) x[2]))
#   
#   ## exon stats
#   str = as.character(strand(exongr)[1])
#   
#   if(str == "+"){
#     start = start(exongr)[1]
#     coordOnGene = pos - start + 1
#   }
#   if(str == "-"){
#     start = end(exongr)[length(exongr)]
#     coordOnGene = start - pos + 1
#   }
#   RLD = log2(sum(frac[,2] * coordOnGene)/sum(frac[,1] * coordOnGene))
#   return(RLD)
# }

