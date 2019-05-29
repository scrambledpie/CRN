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

safe_SE = function(rr){rr = rr[!is.na(rr)]; if(length(rr)==1) return(0) else return(2*sd(rr)/sqrt(length(rr)))}

safe_ME = function(rr){rr = rr[!is.na(rr)]; mean(rr)}

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

# given a list of results and array of Xid, gets the mean and var of each Xid
Get_multiple_OC = function(R, M, with_max=F){
  stopifnot(
    length(R)==length(M),
    typeof(R)=="list",
    typeof(M)%in%c("integer", "double")
  )
  
  uM = sort(unique(M))
  
  lapply(uM, function(mi)Get_OC_mean_var(R[M==mi], with_max=with_max))
  
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



# given a set of OC results matrices, plots the whole lot on a single axes.
Res_plot = function(R, cols=2:100, xr=NULL, yr=NULL, xlab="", ylab="", main="", log=""){
  # R: list of N*3 matrices, each with cols N, OC, and se(OC)
  
  if(is.null(xr)) xr = Get_xy_ran(R)$xr
  if(is.null(yr)) yr = Get_xy_ran(R)$yr
  
  Empty_plot(xr, yr, xlab=xlab, ylab=ylab, main=main, log=log)
  
  k=1
  for(Ri in R){Plot_mean_er(Ri, cols[k]); k=k+1}
  
}