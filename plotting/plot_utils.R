LinInterp = function(X, Y){
  
  O = order(X)
  X = X[O]
  Y = Y[O]
  top_X = X[length(X)]
  
  interp = function(x){
    
    if(x==top_X) return(Y[length(X)])
    if(x==X[1]) return(Y[1])
    if(x<X[1]|x>top_X)return(NA)
    
    l = max(which(X<=x))
    u = min(which(x<=X))
    if(l==u){
      return(Y[u])
    }else{
      return( Y[l] + (Y[u]-Y[l])*(x-X[l])/(X[u]-X[l]) )
    }
  }
  
  return(function(XX)sapply(XX, interp))
}

plotpoly = function(X,Y,ER,coli,bottom=-Inf){
  
  O = order(X)
  X = X[O]
  Y = Y[O]
  
  xx = c(X,rev(X))
  yy = c(Y+ER,rev(Y-ER))
  yy = pmax(yy,bottom)
  # yy = pmax(yy,max(0.0001,min(yy)))
  polygon(xx,yy,col=adjustcolor(coli,alpha=0.5),border = NA)
}

safe_SE = function(rr){rr = rr[!is.na(rr)]; if(length(rr)<2) return(0) else return(2*sd(rr)/sqrt(length(rr)))}

safe_ME = function(rr){rr = rr[!is.na(rr)]; mean(rr)}

safe_range = function(rr){rr = rr[!is.na(rr)]; range(rr)}

# get the mean and var OC of a set of results, returns [X, Y, se(Y)] matrix
Get_OC_mean_var = function(R, with_max=F){
  Nran = sapply(R, function(r)range(r$Cost$N))
  Nvals= round(Nran[1]):round(Nran[2])
  if(with_max){
    OC   = sapply(R, function(r)r$maxY_f - LinInterp(r$Cost$N, r$Cost$P)(Nvals))
  }else{
    OC   = sapply(R, function(r)LinInterp(r$Cost$N, r$Cost$P)(Nvals))
  }
  R_mean = apply(OC,1, safe_ME)
  R_ser  = apply(OC,1, safe_SE)
  cbind(Nvals, R_mean, R_ser)
}
Get_multiple_OC = function(R, M, with_max=F){
  stopifnot(
    length(R)==length(M),
    typeof(R)=="list",
    typeof(M)%in%c("integer", "double")
  )
  
  uM = sort(unique(M))
  
  lapply(uM, function(mi)Get_OC_mean_var(R[M==mi], with_max=with_max))
  
}

# seed reuse time series
Get_seed_reuse = function(R){
  
  Nran = range(sapply(R, function(r)range(r$Cost$N)))
  N0 = Nran[1]
  sc = ncol(R[[1]]$x)
  
  # convert each time series to binary
  bin_stream=function(S){
    top = Nran[2]
    S = S[1:min(length(S), top)]
    out = sapply(Nran[1]:length(S), function(i)S[i]%in%S[1:(i-1)])
    c(out, rep(NA, top-length(S)))
  }
  
  bin_streams = sapply(R, function(R)bin_stream(R$x[,sc]))
  
  ME = apply(bin_streams, 1, safe_ME)
  SE = apply(bin_streams, 1, safe_SE)
  
  return(cbind(Nran[1]:Nran[2], ME, SE))
}
Get_multiple_seed_reuse = function(R, M){
  stopifnot(
    length(R)==length(M),
    typeof(R)=="list",
    typeof(M)%in%c("integer", "double")
  )
  
  uM = sort(unique(M))
  
  lapply(uM, function(mi)Get_seed_reuse(R[M==mi]))
}

# barchart of seed frequencies
Get_seed_bars = function(R, topseed=500){
  
  N0 = min(sapply(R, function(Ri)min(Ri$Cost$N)))

  seeds = lapply(R, function(Ri)Ri$x[1:N0,ncol(Ri$x)])
  get_freqs = function(Seeds)sapply(1:topseed, function(si)sum(si==Seeds))
  Freqs = sapply(seeds, get_freqs)
  mean_F0 = apply(Freqs, 1, mean)
  
  seeds = lapply(R, function(Ri)Ri$x[-(1:N0),ncol(Ri$x)])
  get_freqs = function(Seeds)sapply(1:topseed, function(si)sum(si==Seeds))
  Freqs = sapply(seeds, get_freqs)
  mean_F1 = apply(Freqs, 1, mean)
  
  
  mean_F = rbind(mean_F0, mean_F1)
  # non_empty  = apply(mean_F, 2, sum)>0 
  # mean_F = mean_F[, non_empty]
  
  return(mean_F)
  
}
Get_seed_pairs_singles = function(R, topseed=500){
  
  N0 = min(sapply(R, function(Ri)min(Ri$Cost$N)))
  
  seeds = lapply(R, function(Ri)Ri$x[-(1:N0),ncol(Ri$x)])
  get_freqs = function(Seeds)sapply(1:topseed, function(si)sum(si==Seeds))
  Freqs = sapply(seeds, get_freqs)
  
  Freqs_to_pairs = function(freq)c(sum(freq==1), 2*sum(freq==2))
  singles_pairs = apply(Freqs, 2, Freqs_to_pairs)
  
  
  mean_PW = apply(singles_pairs, 1, mean)
  # non_empty  = apply(mean_F, 2, sum)>0 
  # mean_F = mean_F[, non_empty]
  
  return(mean_PW)
}

# Get the ranges of matrices
Get_xy_ran = function(ResList){
  # ResList is a list of matrices from Get_OC_mean_var
  xran = range(sapply(ResList, function(r)range(r[,1])))
  yran = range(sapply(ResList, function(r)range(r[,2])))
  return(list(xr=xran, yr=yran))
}

# These fucntions initialize a plot and append curves to it
Empty_plot=function(xran, yran, xlab="", ylab="", main="", log=""){
  plot(xran, yran, xlab=xlab, ylab=ylab, col="white", main=main, log=log)
}
Plot_mean_er = function(Res_matrix, col=2){
  lines(Res_matrix[,1], Res_matrix[,2], col=col)
  plotpoly(Res_matrix[,1], Res_matrix[,2], Res_matrix[,3], col=col)
}
Plot_pnts_bars = function(Res_matrix, col=2){
  x = Res_matrix[,1]
  avg = Res_matrix[,2]
  sdev = Res_matrix[,3]
  points(x, avg, col=col, pch=19)
  lines(x, avg, col=col)
  arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3, col=col)
}

plot_errors = function(R, M, P, N=NULL, cols = 2:100){
  
  # R list of matrices
  # M vector of experiment ID
  # P vector of parameter values to plot over
  # N the row in the matrix to pick from
  # cols
  
  stopifnot(
    typeof(R)=="list",
    typeof(M)%in%c("double", "integer"),
    typeof(P)%in%c("double", "integer"),
    length(R)==length(P),
    length(P)==length(M),
    N <= max(sapply(R, nrow)),
    length(unique(M))<=length(cols)
  )
  
  uM = sort(unique(M))
  uP = sort(unique(P))
  
  get_output=function(ri, N) tryCatch(ri[N,2:3], error=function(e)c(NA,NA))
  
  ResList = lapply(uM, function(mi)t(sapply(uP, function(pi) get_output(R[M==mi&P==pi][[1]], N))))
  
  yran = range(sapply(ResList, function(r)safe_range(c(r[,1]+r[,2], r[,1]-r[,2]))))
  
  print(yran)
  
  Empty_plot(range(P), yran)
  
  
  for(i in 1:length(uM)){
    Plot_pnts_bars( cbind(uP, ResList[[i]]) , col=cols[i])
  }
  
  
}

Seeds_barplot = function(R1=Results){
  
  M = sapply(R1, function(r)r$method)
  
  sc = ncol(R1[[1]]$x)
  N0 = min(R1[[1]]$Cost$N)
  
  
  if(sum(M==5)>0){
    get_pairs = function(r){
      SS = r$x[-(1:N0), sc]
      num_seeds = length(unique(SS))
      num_pairs = 2*(length(SS) - num_seeds)
      num_singles = num_seeds - num_pairs
      c(num_singles, num_pairs)
    }
    PWR = R1[M==5]
    SS = apply(sapply(PWR, get_pairs),1,mean)
    cat(SS)
    barplot(SS, names.arg = c("singles", "pairs"), main="CompSph")
    
  }
  
  if(sum(M==7)>0){
    get_pairs = function(r){
      SS = r$x[-(1:N0), sc]
      num_seeds = length(unique(SS))
      num_pairs = 2*(length(SS) - num_seeds)
      num_singles = num_seeds - num_pairs
      c(num_singles, num_pairs)
    }
    PWR = R1[M==7]
    SS = apply(sapply(PWR, get_pairs),1,mean)
    cat(SS)
    barplot(SS, names.arg = c("singles", "pairs"), main="CS+Wiggles")
    
  }
  
  
  
  if(sum(M==4)>0){
    get_seed_freqs = function(r){
      SS = r$x[-(1:N0), sc]
      sapply(1:500, function(s)sum(SS==s))
    }
    get_seed_freqs_start = function(r){
      SS = r$x[1:N0, sc]
      sapply(1:100, function(s)sum(SS==s))
    }
    CRNR = R1[M==4]
    SS = apply(sapply(CRNR, get_seed_freqs), 1, mean)
    top = max(which(SS>0))+1
    SS = SS[1:top]
    
    SSI = apply(sapply(CRNR, get_seed_freqs_start), 1, mean)
    SSI = SSI[1:top]
    
    barplot(rbind(SSI,SS-SSI) , names.arg=paste(1:length(SS)), main = "CompSph")
    
  }
  
  if(sum(M==6)>0){
    get_seed_freqs = function(r){
      SS = r$x[-(1:N0), sc]
      sapply(1:500, function(s)sum(SS==s))
    }
    get_seed_freqs_start = function(r){
      SS = r$x[1:N0, sc]
      sapply(1:500, function(s)sum(SS==s))
    }
    CRNR = R1[M==6]
    SS = apply(sapply(CRNR, get_seed_freqs), 1, mean)
    top = max(which(SS>0))+1
    SS = SS[1:top]
    
    SSI = apply(sapply(CRNR, get_seed_freqs_start), 1, mean)
    SSI = SSI[1:top]
    
    barplot(rbind(SSI,SS-SSI) , names.arg=paste(1:length(SS)), main="CS+wiggles")
    
  }
}


# given a set of OC results matrices, plots the whole lot on a single axes.
Res_plot = function(R, cols=2:100, xr=NULL, yr=NULL, xlab="", ylab="", main="", log=""){
  # R: list of N*3 matrices, each with cols N, OC, and se(OC)
  
  if(is.null(xr)) xr = Get_xy_ran(R)$xr
  if(is.null(yr)) yr = Get_xy_ran(R)$yr
  
  Empty_plot(xr, yr, xlab=xlab, ylab=ylab, main=main, log=log)
  
  k=1
  for(Ri in R){Plot_mean_er(Ri, cols[k]); k=k+1}
  
}





