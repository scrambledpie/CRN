# rm(list=ls())
.libPaths(c(.libPaths(),"~/R/"))
library(MASS)
library(R6)
library(FastGP)
library(mcmc)

source('CRN_BO/utils.R')

Rcpp::sourceCpp('CRN_BO/kernel.cpp')

cat(" kernel.cpp compilation complete. ")


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
##  MODEL

{
  MU    = function(HP, x, xd, yd){
    
    dims = ncol(xd)
    kernel = function(x1,x2){
      SQ = lapply(1:dims,function(d)outer(x1[,d],x2[,d],"-")^2)
      
      HP[dims+1]*exp(-0.5*Reduce("+",mapply("*",SQ,1/HP[1:dims]^2,SIMPLIFY = F)))
    }
    mean(yd)+kernel(x,xd)%*%rcppeigen_invert_matrix(kernel(xd,xd)+HP[dims+2]*diag(length(yd)))%*%(yd-mean(yd))
    
  }
  LhoodBuilder  = function(xd, yd, with_prior=F, XRAN=NULL){
    
    stopifnot(
      length(xd)%%length(yd)==0,
      all(!is.na(xd)),
      all(!is.na(yd)),
      all(is.finite(xd)),
      all(is.finite(yd))
    )
    
    if(with_prior){
      cat(" Lhoodbuilding with prior, ")
      LX_PriorMean = XRAN[2,] - XRAN[1,]
      LX_PriorPrec = 0.25/(LX_PriorMean)^2
      
      SX_PriorMean = var(yd)
      SX_PriorPrec = 0.25/(SX_PriorMean)^2
      
      NV_PriorRate = 1/(1000000 * SX_PriorMean * 0.5)
      
      dims = ncol(xd)
      
      H_prior = function(HP){
        
        # length scale prior log density
        PLX = sum(-0.5*(HP[1:dims]-LX_PriorMean)^2*LX_PriorPrec)
        
        # signal variance prior log density
        PSX = -0.5*(HP[dims+1]-SX_PriorMean)^2*SX_PriorPrec
        
        # noise parameter prior density
        PNV = -NV_PriorRate * sum(HP[dims+2])
        
        PLX + PSX + PNV
        # browser()
      }
      DH_prior = function(HP){
        
        # signal variance prior log density
        DLX = -(HP[1:dims]-LX_PriorMean)*LX_PriorPrec
        
        # signal variance prior log density
        DSX = -(HP[dims+1]-SX_PriorMean)*SX_PriorPrec
        
        # noise parameter prior density
        DNV = -NV_PriorRate
        
        c(DLX, DSX, DNV) * HP
      }
      
    }else{
      cat(" Lhoodbuilding without prior, ")
      H_prior = function(HP)0
      DH_prior = function(HP)0
    }
    
    if(all(xd[,ncol(xd)]%%1==1) )dims = ncol(xd)-1 else dims = ncol(xd)
    GPmean0 <<- mean(yd)
    
    SQ = lapply(1:dims,function(d)outer(xd[,d],xd[,d],"-")^2)
    ydm = yd - GPmean0
    NX = nrow(xd)
    alpha = 0
    KK = 0
    SE = 0
    lastLTheta = rep(0, length(SQ)+2)
    History = matrix(0,0,length(SQ)+3)
    lhood_term3 = 0.5*length(yd)*log(2*pi)
    
    
    # define likelihood and its gradient
    
    Obj = function(LTheta){
      
      
      if(length(LTheta)!=length(SQ)+2) stop("wrong number of hpars, got ", length(LTheta), " expected ", length(SQ)+2,"\n")
      
      lastLTheta <<- LTheta
      Theta= exp(LTheta)
      iL = -0.5/(Theta[1:dims]^2)
      SE <<-exp(Reduce("+",mapply("*",SQ,iL,SIMPLIFY = F)))
      KK <<- Theta[dims+1]*SE + Theta[dims+2]*diag(NX)
      
      minl = Inf
      
      inv_KK <<- rcppeigen_invert_matrix(KK)
      alpha <<- ydm %*% inv_KK
      diagK <- rcppeigen_get_diag(KK)
      out =  - sum(log(diagK)) - (1/2) * sum(alpha  * ydm) - lhood_term3 + H_prior(Theta)
      
      # cat("prior density is ",H_prior(Theta), " " )
      if(is.na(out)) return(minl) else if(out<minl)minl<<-out
      
      History<<- rbind(History,c(Theta,out))
      # browser()
      return(out[1])
    }
    
    DObj = function(LTheta){
      
      if(length(LTheta)!=length(SQ)+2) stop("wrong number of hpars, got ", length(LTheta), " expected ", length(SQ)+2,"\n")
      
      if(any(LTheta!=lastLTheta)) dead= Obj(LTheta)
      
      # cat("size of alpha", length(alpha))
      # cat("size of invKK", dim(inv_KK))
      
      Theta = exp(LTheta)
      il3 = 1/(Theta[1:dims]^3)
      
      alpha = as.numeric(alpha)
      D0 = outer(alpha,alpha) - inv_KK
      # cat(length(alpha), dim(inv_KK))
      
      Dl1 =sapply(1:dims,function(d) sum(D0 * SQ[[d]]*SE)*il3[d]*Theta[dims+1])
      DS2 = sum(D0 * SE)
      DV2 = sum(diag(D0))
      
      c(Dl1,DS2,DV2)*Theta*0.5 + DH_prior(Theta)
      
    }
    
    if (length(yd)==1){
      top=NULL
      low=NULL
    }else{
      xtop = 1.5*(apply(xd[,1:dims,drop=F],2,max) - apply(xd[,1:dims,drop=F],2,min))
      ytop = 0.5*(max(yd)-min(yd))^2
      log_top = log( c(xtop, ytop, ytop) )
    
      xlow = xtop/100
      ylow = ytop/1000
      log_low  = log( c(xlow, ylow, ylow/100))
    }
    
    bestpars = function(){
      History
    }
    
    return(list(LH=Obj, DLH=DObj, ParsHist=bestpars, lo=log_low, hi=log_top))
  }
  LhoodBuilder_standard  = function(xd, yd, with_prior=F, XRAN=NULL){
    FUNS = LhoodBuilder(xd, yd, with_prior, XRAN)
    
    FLH = FUNS$LH
    FDLH = FUNS$DLH
    
    Obj = function(Hpars){
      # print(Hpars)
      FLH(log(Hpars))
    }
    
    DObj = function(Hpars){
      FDLH(log(Hpars))/Hpars
    }
    
    low = exp(FUNS$lo)
    top = exp(FUNS$hi)
    
    return(list(LH=Obj, DLH=DObj, ParsHist=FUNS$ParsHist, lo=low, hi=top))
  }
  CRNLHood = R6Class("LhoodOptimizer",public = list(
    npars  = NULL,
    xd     = NULL,
    yd     = NULL,
    yd_o   = NULL,
    ydm    = NULL,
    ymean  = NULL,
    alpha  = NULL,
    SQ     = NULL,
    Seeds  = NULL,
    dims   = NULL,
    SEt    = NULL,
    SEe    = NULL,
    Noise  = 0.1,
    KK     = NULL,
    invK   = NULL,
    lastLHpars = NULL,
    lastHpars = NULL,
    EvalHistory=NULL,
    Hpar_History = list(),
    htop   = NULL,
    hlow   = NULL,
    ltop   = NULL,
    llow   = NULL,
    HP     = NULL,
    XRAN   = NULL,
    iKY    = NULL,
    iK     = NULL,
    LX_PriorMean = NULL,
    LX_PriorPrec = NULL,
    SX_PriorMean = NULL,
    SX_PriorPrec = NULL,
    NV_PriorRate = NULL,
    lhood_term3  = NULL,
    copula       = FALSE,
    
    initialize = function(xd, yd, XRAN, copula=F){
      
      self$copula  = copula
      self$dims    = ncol(xd)-1
      self$npars   = self$dims*2 + 3
      self$xd      = xd
      self$yd_o    = yd
      
      if(self$copula){
        self$yd      = qnorm(order(yd)/(length(yd)+1))
      }else{
        self$yd      = self$yd_o
      }
      
      self$XRAN    = XRAN

    },
    
    Lhood_Prep = function(){
      cat("prepping Lhood, ")
      
      self$LX_PriorMean = XRAN[2,] - XRAN[1,]
      self$LX_PriorPrec = 0.25/(self$LX_PriorMean)^2
      
      self$SX_PriorMean = var(self$yd)
      self$SX_PriorPrec = 0.25/(self$SX_PriorMean)^2
      
      self$NV_PriorRate = 1/(1000000 * self$SX_PriorMean * 0.5)
      
      if(!is.null(self$HP)){
        self$lastLHpars = log(self$HP)
        self$lastHpars  = self$HP
      }
      self$SQ    = lapply(1:self$dims, function(d) -0.5*outer(self$xd[,d], self$xd[,d], "-")^2)
      self$Seeds = outer(self$xd[,ncol(self$xd)], self$xd[,ncol(self$xd)], "==")
      self$ymean = mean(self$yd)
      self$ydm   = self$yd - mean(self$yd)
      
      xtop = 1.5*(self$XRAN[2,]-self$XRAN[1,])
      ytop = 0.1*(max(self$yd)-min(self$yd))^2
      self$htop = c(xtop, ytop, xtop, ytop, ytop, ytop)
      self$ltop = log(self$htop)
      
      xlow = xtop/200
      ylow = ytop/1000
      self$hlow = c(xlow, ylow, xlow, ylow/10, ylow/100, ylow/100)
      self$llow = log(self$hlow)
      
      self$EvalHistory = NULL
      self$Noise = 0.000*diag(length(self$yd))
      self$lhood_term3 = 0.5*length(self$yd)*log(2*pi)
    },
    
    LHyperPrior = function(Hpars){
      
      dd = self$dims 
      if(length(Hpars)!=dd*2+4)stop("LHprior wrong input length hpars, get ", length(Hpars))
      # length scale prior log density
      PLX = sum(-0.5*(Hpars[1:dd]-self$LX_PriorMean)^2*self$LX_PriorPrec)
      
      # wiggle length scale prior log density
      PLXw = sum(-0.5*(Hpars[1+dd+1:dd]-self$LX_PriorMean*0.5)^2*self$LX_PriorPrec)
      
      
      # signal variance prior log density
      PSX = -0.5*(Hpars[dd+1]-self$SX_PriorMean)^2*self$SX_PriorPrec
      
      # wiggle variance prior log density
      PSXw = -0.5*(Hpars[2*(dd+1)]-self$SX_PriorMean*0.5)^2*self$SX_PriorPrec
      
      # offset variance prior log density
      PC = -0.5*(Hpars[2*(dd+1)+1]-self$SX_PriorMean*0.5)^2*self$SX_PriorPrec
      
      # white noise parameter prior density
      PNV = -self$NV_PriorRate * Hpars[length(Hpars)]
      
      return( PLX + PLXw + PSX + PSXw + PC + PNV )
  
    },
    DLHyperPrior = function(Hpars){
      
      dd = self$dims
      
      if(length(Hpars)!=dd*2+4)stop("DLHprior wrong input length hpars, get ", length(Hpars))
      # grad length scales density
      DLX = -(Hpars[1:dd]-self$LX_PriorMean)*self$LX_PriorPrec

      # grad wiggle length scales density
      DLXw = -(Hpars[1+dd+1:dd]-self$LX_PriorMean*0.5)*self$LX_PriorPrec
      
      # grad signal variance density
      DSX = -(Hpars[dd+1]-self$SX_PriorMean)*self$SX_PriorPrec
      
      # grad wiggle signal variance density
      DSXw = -(Hpars[2*(dd+1)]-self$SX_PriorMean*0.5)*self$SX_PriorPrec
      
      DC  =  -(Hpars[2*(dd+1)+1]-self$SX_PriorMean*0.5)*self$SX_PriorPrec
      
      # grad signal variance density
      DNV = -self$NV_PriorRate
      
      c(DLX, DSX, DLXw, DSXw, DC, DNV)
      
    },
    
    Lhood  = function(LHpars){
      
      dd = self$dims
      
      if (is.null(self$SQ)) self$Lhood_Prep()
      if (length(self$yd)!=nrow(self$SQ[[1]])) self$LhoodPrep()
      
      if (length(LHpars)!=self$dims*2+4) stop("lhood wrong input length hpars, get ", length(LHpars))
      
      self$lastLHpars = LHpars
      Hpars = exp(LHpars)
      Di   = 1:dd
      LXt  = Hpars[Di]
      S2t  = Hpars[dd+1]
      LXe  = Hpars[dd+1 + Di]
      S2e  = Hpars[dd+1 + dd+1]
      C2e  = Hpars[dd+1 + dd+1 + 1]  
      V2e  = Hpars[length(Hpars)]
      
      
      SEt = exp(Reduce("+",mapply("*",self$SQ,(1/LXt^2),SIMPLIFY = F)))
      SEe = exp(Reduce("+",mapply("*",self$SQ,(1/LXe^2),SIMPLIFY = F)))
      
      # browser()
      
      KK    = S2t*SEt + (S2e*SEe + C2e)*self$Seeds + V2e*diag(length(self$yd))
      invK  = rcppeigen_invert_matrix(KK)
      
      if(any(is.nan(invK)) ) stop("Nans in matrix")
      
      self$SEt  = SEt
      self$SEe  = SEe
      self$KK   = KK
      self$invK = invK
      
      self$alpha = as.numeric(self$ydm%*%self$invK)
      
      OO = -sum(log(rcppeigen_get_diag(self$KK))) -0.5*sum(self$ydm*self$alpha) - self$lhood_term3
      
      
      if(all(is.finite(c(Hpars,OO)))) {
        self$EvalHistory = rbind(self$EvalHistory, c(Hpars,OO))
      }
      
      OO
      
    },
    DLhood = function(LHpars){
      dd = self$dims
      
      if(length(LHpars)!=dd*2+4)stop("Dlhood wrong input length hpars, get ", length(LHpars))
      
      if (is.null(self$SQ)) self$Lhood_Prep()
      if (length(self$yd)!=nrow(self$SQ[[1]])) self$LhoodPrep()
      
      if(any(LHpars!=self$lastLHpars)) BB = self$Lhood(LHpars)
      
      Di = 1:self$dims
      Hpars = exp(LHpars)
      
      LXt  = Hpars[Di]
      S2t  = Hpars[dd+1]
      LXe  = Hpars[dd+1 + Di]
      S2e  = Hpars[dd+1 + dd+1]
      C2e  = Hpars[dd+1 + dd+1 +1]
      
      DD   = outer(self$alpha,self$alpha) - self$invK
      DLXt = -2*sapply(Di,function(d) sum(DD*self$SQ[[d]]*self$SEt)*(S2t/LXt[d]^3))
      DS2t = sum(DD*self$SEt)
      DLXe = -2*sapply(Di,function(d) sum(DD*self$SQ[[d]]*self$SEe*self$Seeds)*(S2e/LXe[d]^3))
      DS2e = sum(DD*self$SEe*self$Seeds)
      DC2e = sum(DD*self$Seeds)
      DV2e = sum(diag(DD))
      
      return(0.5*c(DLXt,DS2t,DLXe,DS2e,DC2e,DV2e)*Hpars)
      
    },
    Lhood_standard  = function(Hpars){
      
      dd = self$dims
      
      if (is.null(self$SQ)) self$Lhood_Prep()
      if (length(self$yd)!=nrow(self$SQ[[1]])) self$Lhood_Prep()
      
      if (length(Hpars)!=self$dims*2+4) stop("lhood wrong input length hpars, get ", length(Hpars))
      
      self$lastHpars = Hpars
      Di   = 1:dd
      LXt  = Hpars[Di]
      S2t  = Hpars[dd+1]
      LXe  = Hpars[dd+1 + Di]
      S2e  = Hpars[dd+1 + dd+1]
      C2e  = Hpars[dd+1 + dd+1 + 1]  
      V2e  = Hpars[length(Hpars)]
      
      
      SEt = exp(Reduce("+",mapply("*",self$SQ,(1/LXt^2),SIMPLIFY = F)))
      SEe = exp(Reduce("+",mapply("*",self$SQ,(1/LXe^2),SIMPLIFY = F)))
      
      # browser()
      
      KK    = S2t*SEt + (S2e*SEe + C2e)*self$Seeds + V2e*diag(length(self$yd))
      invK  = rcppeigen_invert_matrix(KK)
      
      if(any(is.nan(invK)) ) {cat("K is nan");stop("Nans in matrix")}
      
      self$SEt  = SEt
      self$SEe  = SEe
      self$KK   = KK
      self$invK = invK
      self$alpha = as.numeric(self$ydm%*%self$invK)
      
      OO = -sum(log(rcppeigen_get_diag(self$KK))) - 0.5*sum(self$ydm*self$alpha) - self$lhood_term3
      
      if(all(is.finite(c(Hpars,OO)))) {
        self$EvalHistory = rbind(self$EvalHistory, c(Hpars,OO))
      }
      
      OO
      
    },
    DLhood_standard = function(Hpars){
      
      dd = self$dims
      
      if(length(Hpars)!=dd*2+4)stop("Dlhood wrong input length hpars, get ", length(Hpars))
      
      if (is.null(self$SQ)) self$Lhood_Prep()
      if (length(self$yd)!=nrow(self$SQ[[1]])) self$Lhood_Prep()
      
      if(any(Hpars!=self$lastHpars)) BB = self$Lhood_standard(Hpars)
      
      Di = 1:self$dims
      
      LXt  = Hpars[Di]
      S2t  = Hpars[dd+1]
      LXe  = Hpars[dd+1 + Di]
      S2e  = Hpars[dd+1 + dd+1]
      C2e  = Hpars[dd+1 + dd+1 +1]
      V2e  = Hpars[length(Hpars)]
      
      DD   = outer(self$alpha,self$alpha) - self$invK
      DLXt = -2*sapply(Di,function(d) sum(DD*self$SQ[[d]]*self$SEt)*(S2t/LXt[d]^3))
      DS2t = sum(DD*self$SEt)
      DLXe = -2*sapply(Di,function(d) sum(DD*self$SQ[[d]]*self$SEe*self$Seeds)*(S2e/LXe[d]^3))
      DS2e = sum(DD*self$SEe*self$Seeds)
      DC2e = sum(DD*self$Seeds)
      DV2e = sum(diag(DD))
      
      0.5*c(DLXt,DS2t,DLXe,DS2e,DC2e,DV2e)
      
    },
    
    BitbyBitOptimise = function(){
      
      Dead = self$LhoodPrep()
      
      if(!is.null(self$lastLHpars)){
        cat("Sinoptim prev pars, ")
        dead = exp(Optimizer(self$Lhood,
                             self$DLhood,
                             rbind(self$llow, self$ltop),
                             x0=self$lastLHpars, 
                             Ns=0 ))
      }
      
      # make some heuristic bounds on hpars
      xtop = 1.5*(apply(self$xd[,-ncol(self$xd),drop=F],2,max) - apply(self$xd[, -ncol(self$xd), drop=F],2,min))
      ytop = 0.5*(max(self$yd)-min(self$yd))^2
      top = c(xtop,ytop,xtop,ytop,ytop)
      self$htop = top
      self$ltop = log(top)
      
      xlow = xtop/100
      ylow = ytop/1000
      low  = c(xlow, ylow, xlow, ylow/10, ylow/10) 
      self$hlow = low
      self$llow = log(low)
      
      
      XXi   = self$xd
      X     = self$xd[,1:self$dims,drop=F]
      Xt    = t(X)
      
      
      ####  estimate seed specific Hpars [[ Se, Sc ]]
      Seeds = self$xd[,self$dims+1]
      
      # estimate varaince of offsets
      Sc    = var( sapply(unique(Seeds),function(s)mean(self$yd[Seeds==s])))
      
      # gather x neighbours of different seeds
      NN    = apply(X,1,function(x)apply((x-Xt)^2,2,sum)<0.01) & ((outer(Seeds,Seeds,"!=")| diag(length(Seeds))==1))
      
      # be sure to use each y value only once, first exclude those with no neighbours, then get average var of each cluster
      unused_y = apply(NN,1,sum)>1
      Se0   = 0
      k     = 0
      for( i in 1:nrow(NN)){ if(unused_y[i]){
        k=k+1
        Se0 = Se0 + var(self$yd[NN[i,]])
        unused_y[NN[i,]]=F
      }}
      
      # browser()
      
      # estimate wiggle signal variance
      Se0   = abs(Se0/k - Sc)
      
      #  centre y values from each seed, fit standard GP to all data with "known" wiggle noise variance to get [[ Lt , St ]]
      Ym         = self$yd
      for(s in unique(Seeds)) Ym[Seeds==s] = Ym[Seeds==s] - mean(Ym[Seeds==s])
      Ltf        = LhoodBuilder(X,Ym)

      cat("theta, ")
      LStheta    = Freezer(Ltf$LH, Ltf$DLH, rbind(Ltf$lo, Ltf$hi), self$dims+2, log(Se0))$xmax
      
      
      # get residuals between theta-y, fit indep GPs to each seed wiggle to get [[ Le , Se ]]
      MM        = MU(c(LStheta, Se0), X, X, Ym)
      Ym2       = Ym-MM
      Lhoodfuns = lapply(unique(Seeds),function(s)LhoodBuilder(X[Seeds==s,,drop=F],Ym2[Seeds==s]))
      LesLhood  = function(lt)sum(sapply(Lhoodfuns,function(L)L$LH(lt)))
      LesDLhood = function(lt){DLH = sapply(Lhoodfuns,function(L)L$DLH(lt)); apply(DLH,1,sum)}
      
      lowest  = Lhoodfuns[[1]]$lo
      highest = Lhoodfuns[[1]]$hi
      if(length(Lhoodfuns)>1){
        for (i in 2:length(Lhoodfuns)){
          lowest    = cbind(lowest,  Lhoodfuns[[i]]$lo)
          highest   = cbind(highest, Lhoodfuns[[i]]$hi)
        }
        lowest = apply(lowest, 1, min)
        highest = apply(highest, 1, max)
      }
      
      cat("epsilon (assuming noisy GPs), ")
      LSe       = Optimizer(LesLhood,
                            LesDLhood,
                            rbind(lowest, highest), 
                            debugging=T)$xmax[1:(self$dims+1)]
      cat("all, ")
      
      # We constrain the sum of variance terms to be half the range of y values.
      TotalVar = (max(self$yd) - min(self$yd))^2/4
      log_OffsetVar = function(log_St, log_Se)log(max(0.01, TotalVar - exp(log_St) - exp(log_Se)))
      dlog_OffsetVar = function(log_St, log_Se){St=exp(log_St); Se=exp(log_Se); -c(St, Se)/max(0.01, TotalVar - St-Se)}
      mod_Lhood = function(LHpars){
        log_St = LHpars[self$dims+1]
        log_Se = LHpars[2*(self$dims+1)]
        allpars = c(LHpars, log_OffsetVar(log_St, log_Se))
        self$Lhood(allpars)
      }
      mod_DLhood = function(LHpars){
        log_St = LHpars[self$dims+1]
        log_Se = LHpars[2*(self$dims+1)]
        allpars = c(LHpars, log_OffsetVar(log_St, log_Se))
        grad = self$DLhood(allpars)
        
        grad_off = dlog_OffsetVar(log_St, log_Se)
        mod_grad = c(grad[1:self$dims],
                     grad[self$dims+1] + grad[length(grad)]*grad_off[1],
                     grad[self$dims+1 + 1:self$dims],
                     grad[2*(self$dims+1)] + grad[length(grad)]*grad_off[2]
                    )
        mod_grad
      }
      mod_low = self$llow[-length(self$llow)]
      mod_top = self$ltop[-length(self$ltop)]
      # browser()
      # finally optimise all parameters together
      # dead = Optimizer(self$Lhood,
      #                  self$DLhood,
      #                  rbind(self$llow, self$ltop),
      #                  x0=c(LStheta, LSe, log(Sc)),
      #                  Ns=0)
      
      # browser()
      dead = Optimizer(mod_Lhood,
                       mod_DLhood,
                       rbind(mod_low, mod_top),
                       x0=c(LStheta, LSe),
                       Ns=0)
      
      
      topHP = self$Best_HP_History()
      
      cat("completed: ", signif(topHP,3), "\n")

      topHP
      
    },
    
    Unique_seeds_Optimise   = function(with_prior=T){

      
      LFUNS = LhoodBuilder(self$xd[,1:self$dims], self$yd, with_prior=with_prior, XRAN = self$XRAN)

      if(!is.null(self$HP)){
        oldHP = c(self$HP[1:(self$dims+1)], sum( self$HP[-2:0+length(self$HP)] ) )
        oldHP = log(oldHP)
      }else{
        oldHP = NULL
      }
      
      topHP = Optimizer_2(LFUNS$LH, LFUNS$DLH, rbind(LFUNS$lo, LFUNS$hi), x0=oldHP)$xmax
      
      topHP = exp(topHP)
      
      topHP = c(topHP[1:(self$dims+1)], rep(0.00001, self$dims), 0, 0, topHP[length(topHP)])
      
      topHP
      
    },
    Unique_seeds_Optimise_finetune   = function(with_prior=T){
      
      if(is.null(self$HP)) stop("cannot fine tune without old Hpars!")
      
      oldHP = c(self$HP[1:(self$dims+1)], sum(self$HP[-2:0+length(self$HP)]))
      oldHP = log(oldHP)
      
      LFUNS = LhoodBuilder(self$xd[,1:self$dims], 
                           self$yd, 
                           with_prior=with_prior, 
                           XRAN = self$XRAN)
      
      topHP = Optimizer_2(LFUNS$LH, 
                          LFUNS$DLH, 
                          rbind(LFUNS$lo, LFUNS$hi), 
                          x0=oldHP, 
                          N0=0,
                          maxevals = 20)$xmax
      
      topHP = exp(topHP)
      
      topHP = c(topHP[1:(self$dims+1)], rep(0.00001, self$dims), 0,0, topHP[length(topHP)])
      
      topHP
    },
    
    HyperOptimise_CompSphere = function(with_prior=T){
      dd = self$dims
      Di = 1:self$dims
      
      # hpar mapping functions
      CS_to_full = function(hp) c(hp[Di], hp[dd+1], rep(0.0001, dd), 0, hp[dd+2:3])
      full_to_CS = function(hp) c(hp[Di], hp[dd+1], hp[2*dd+3:4])
      
      # First optimize a standard GP with i.i.d. noise
      dead  = self$Lhood_Prep()
      topHP = self$Unique_seeds_Optimise(with_prior = with_prior)
      bestV = topHP[length(topHP)]
      cat("unique seed model optimized, ")
      
      # Next make a functions with one set of length scales for both GPs
      topLH  = -Inf
      topHP2  = 0
      CS_lhood = function(hp){
        
        # make full hpar vector
        hp = CS_to_full(hp)
        
        # try to evaluate, otherwise give -10000000
        tryCatch({
          O = self$Lhood_standard(hp) + with_prior*self$LHyperPrior(hp)
        }, error=function(e)-100000)
        
        # store best result
        if(O>topLH){
          topLH  <<- O
          topHP2 <<- hp
        }
        
        return(O)
      }
      CS_grad_lhood = function(hp){
        
        # make full hpar vector
        hp = CS_to_full(hp)
        
        # evauate hpar grads
        out = self$DLhood_standard(hp)+ with_prior*self$DLHyperPrior(hp)
        
        # cherry pick necessary grads
        return(full_to_CS(out))
      }
      
      # if old parrs exists, evaluate them
      if(!is.null(self$HP)){
        hp0 = full_to_CS(self$HP)
        dead = CS_lhood(hp0)
      }
      
      # Optimize the mixture of noise models
      topObj = -Inf
      topRho = 0
      Obj = function(rho){
        if(any(rho<0) | any(rho>1)) return(-1000000)
        
        x = c(topHP[1:(self$dims+1)],
              rho*bestV,
              (1-rho)*bestV)
        
        O = CS_lhood(x)
        
        if(O>topObj){
          topObj <<- O
          topRho <<- rho
        }
        
        return( O )
      }
      dead = optimise(f=Obj, interval = c(0.0001, 0.9999), maximum = T)
      cat("CS Noise optimized,...")
      
      
      # now optimize everything!
      uu   = full_to_CS(self$htop)
      ll   = full_to_CS(self$hlow)
      hp0  = full_to_CS(topHP2)
      dead = Optimizer_2(CS_lhood,
                         CS_grad_lhood,
                         N0    = 0,
                         x0    = hp0,
                         ran   = rbind(ll, uu),
                         maxevals = 100
                        )
      
      cat("Comp Sph Hyperparams optimized! ")
      return(topHP2)
      
    },
    HyperOptimise_CompSphere_finetune = function(with_prior=T){
      dd = self$dims
      Di = 1:dd
      
      CS_to_full = function(hp) c(hp[Di], hp[dd+1], rep(0.0001, dd), 0, hp[dd+2:3])
      full_to_CS = function(hp) c(hp[Di], hp[dd+1], hp[2*dd+3:4])
      
      
      # First optimize a standard GP with i.i.d. noise
      dead  = self$Lhood_Prep()
      
      # Next make lhoods with zero wiggles
      topLH  = -Inf
      topHP2  = 0
      CS_lhood = function(hp){
        
        # make full hpar vector
        hp = CS_to_full(hp)
        
        # try to evaluate, otherwise give -10000000
        tryCatch({
          O = self$Lhood_standard(hp) + with_prior*self$LHyperPrior(hp)
        }, error=function(e)-100000)
        
        # store best result
        if(O>topLH){
          topLH  <<- O
          topHP2 <<- hp
        }
        
        return(O)
      }
      CS_grad_lhood = function(hp){
        
        # make full hpar vector
        hp = CS_to_full(hp)
        
        # evauate hpar grads
        out = self$DLhood_standard(hp)+ with_prior*self$DLHyperPrior(hp)
        
        # cherry pick necessary grads
        return(full_to_CS(out))
        
      }
      
      # if old pars exists, evaluate them
      if(!is.null(self$HP)){
        hp0   = full_to_CS(self$HP)
        dead  = CS_lhood(hp0)
        topHP = self$HP
        bestV = sum(self$HP[length(self$HP) + -2:0])
      }else{stop("CS finetuning: can't finetune without Hpars!")}
      
      # Optimize the mixture of noise models
      topObj = -Inf
      topRho = 0
      Obj = function(rho){
        if(any(rho<0) | any(rho>1)) return(-1000000)
        
        x = c(topHP[1:(dd+1)],
              rho*bestV,
              (1-rho)*bestV)
        
        O = CS_lhood(x)
        
        if(O>topObj){
          topObj <<- O
          topRho <<- rho
        }
        
        return( O )
      }
      dead = optimise(f=Obj, interval = c(0.0001, 0.9999), maximum = T)
      cat("CS Noise optimized,...")
      # browser()
      
      
      # now optimize everything!
      uu   = full_to_CS(self$htop)
      ll   = full_to_CS(self$hlow)
      hp0  = full_to_CS(topHP2)
      dead = Optimizer_2(CS_lhood,
                         CS_grad_lhood,
                         N0    = 0,
                         x0    = hp0,
                         ran   = rbind(ll, uu),
                         maxevals = 100
      )
      cat("Comp Sph Hyperparams optimized! ")
      return(topHP2)
      
    },
    
    HyperWiggle_optimise = function(with_prior=T){
      
      dead  = self$Lhood_Prep()
      U_HP = self$Unique_seeds_Optimise(with_prior = with_prior)
      U_HP = U_HP[c( 1:(self$dims+1), length(U_HP))]
      
      cat("unique seed model optimized, ")
      
      # get the residuals, centre them, fit i.i.d GP
      Res = self$yd - MU(U_HP, self$xd[,1:self$dims], self$xd[,1:self$dims], self$yd)
      
      S = self$xd[,self$dims+1]
      Su = sort(unique(S))
      
      So = sapply(Su, function(si)mean(Res[S==si]))
      Offset_sig = var(So)
      
      W_Lhoods = list()
      Res2 = Res
      t=1
      for (i in 1:length(Su)){
        mask = S==Su[i]
        Res2[mask] = Res[mask] - So[i]
        if(sum(mask)>1){
          xsi = self$xd[mask,1:self$dims]
          yi  = self$yd[mask]
          cat(" xdim",dim(xsi))
          cat(" y",length(yi), "\n")
          W_Lhoods[[t]] = LhoodBuilder_standard(xsi, yi, T, self$XRAN)
          t=t+1
        }
      }
      
      W_Obj = function(HP)sum(sapply(W_Lhoods, function(fi)fi[[1]](HP)))
      W_DObj = function(HP)apply(sapply(W_Lhoods, function(fi)fi[[2]](HP)),1,sum)
      W_hi = apply(sapply(W_Lhoods, function(lh)lh$hi), 1, mean)
      W_lo = apply(sapply(W_Lhoods, function(lh)lh$lo), 1, mean)
      
      
      browser()
      
      
      if(!is.null(self$HP)){
        W_lx_old = self$HP[self$dims + 1 + 1:self$dims]
        W_sig_old = self$HP[2*(1+self$dims)]
      }else{
        W_lx_old = U_HP[1:self$dims]
        W_sig_old = (max(Res2) - min(Res2))^2/4
      }
      
      W_noise = max(U_HP[length(U_HP)] - W_sig_old - Offset_sig, W_sig_old/100)
      
      cat(W_noise)
      
      W_HP0 = c(W_lx_old, W_sig_old, W_noise)
        
      cat(" here we go wiggles!! ")
      W_topHP = Optimizer_2(W_Obj, W_DObj, rbind(W_lo, W_hi), x0=W_HP0)$xmax
      
      cat(" wiggles optimized, ")
      
      # Now we have all the parameters, we can just finetune!
      UW_HP0 = c( U_HP[1:(self$dims+1)], W_topHP[-length(W_HP0)], Offset_sig )
      
      # browser()
      
      topObj  = -Inf
      Total_topHP  = 0
      Total_lhood = function(hp){
        O = self$Lhood_standard(hp) + with_prior*self$LHyperPrior(hp)
        if(O>topObj){
          topObj <<- O
          Total_topHP <<- hp
        }
        return(O)
      }
      
      dead = Total_lhood(UW_HP0)
      
      grad_lhood = function(hp)self$DLhood_standard(hp)+ with_prior*self$DLHyperPrior(hp)
      
      # tryCatch({
      dead = Optimizer_2(Total_lhood,
                         grad_lhood,
                         N0    = 0,
                         x0    = UW_HP0,
                         ran   = rbind(self$hlow, self$htop) 
                         )
      # },finally = function(ee)cat("finetune lhood went bust!"))
      
      cat("completed HPar optimize!")
      
      return(Total_topHP)
      
      
      
      
    },
    HyperWiggle_optimise_finetune = function(with_prior=T){
      
     
      
      # Now we have all the parameters, we can just finetune!
      UW_HP0 = self$HP
      
      topObj  = -Inf
      Total_topHP  = 0
      Total_lhood = function(hp){
        cat("\n Wiggle hood")
        O = self$Lhood_standard(hp) + with_prior*self$LHyperPrior(hp)
        cat(O)
        if(O>topObj){
          topObj <<- O
          Total_topHP <<- hp
        }
        return(O)
      }
      
      dead = Total_lhood(UW_HP0)
      
      grad_lhood = function(hp)self$DLhood_standard(hp)+ with_prior*self$DLHyperPrior(hp)
      
      tryCatch({
      dead = Optimizer_2(Total_lhood,
                         grad_lhood,
                         N0    = 0,
                         x0    = UW_HP0,
                         ran   = rbind(self$hlow, self$htop) 
      )}, error=function(e)cat("\n\n!!!!!!!!!!!!\ntuner busted, stopped early\n!!!!!!!!!!!!!!\n\n"))
      
      
      cat("completed HPar wiggle finetuning")
      
      return(Total_topHP)
      
    },
    
    HyperOptimise_fixedwiggle = function(with_prior=T){
      Di = 1:self$dims
      
      # First optimize a standard GP with i.i.d. noise
      dead  = self$Lhood_Prep()
      topHP = self$Unique_seeds_Optimise(with_prior = with_prior)
      bestV = topHP[length(topHP)]
      cat("unique seed model optimized, ")
      
      # Next make functions with one set of length scales for both GPs
      topLH  = -Inf
      topHP2  = 0
      lhood = function(hp){
        hp = c(hp[Di], hp[self$dims+1], hp[Di], hp[self$dims+2:4])
        tryCatch({
          O = self$Lhood_standard(hp) + with_prior*self$LHyperPrior(hp)
        }, error=function(e)-100000)
        if(O>topLH){
          topLH  <<- O
          topHP2 <<- hp
        }
        return(O)
      }
      grad_lhood = function(hp){
        hp = c(hp[Di], hp[self$dims+1], hp[Di], hp[self$dims+2:4])
        out = self$DLhood_standard(hp)+ with_prior*self$DLHyperPrior(hp)
        out[Di] = out[Di] + out[self$dims+1+Di]
        c(out[Di], out[self$dims+1], out[2*self$dims+2:4])
      }
      if(!is.null(self$HP)){
        hp0 = self$HP
        hp0 = c(hp0[Di], hp0[self$dims+1], hp0[2*self$dims+2:4])
        dead = lhood(hp0)
      }
      
      # Optimize the mixture of noise models
      topObj = -Inf
      topRho = 0
      Obj = function(rho){
        if(any(rho<0) | any(rho>1)) return(-1000000)
        x = c(topHP[1:(self$dims+1)],
              (1-rho[2])*rho[1]    *bestV,
              (1-rho[2])*(1-rho[1])*bestV,
              rho[2]*bestV)
        
        O = lhood(x)
        
        if(O>topObj){
          topObj <<- O
          topRho <<- rho
        }
        
        return( O )
      }
      dead = apply(LHCdim(40,2), 1, Obj)
      dead = optim(f=Obj, par=topRho, method="Nelder-Mead", control = list(maxit=25, fnscale=-1))
      cat("Noise model mixture optimized,...")
      # browser()
      
      
      # now optimize everything!
      uu = self$htop
      uu = c(uu[Di], uu[self$dims+1], uu[2*self$dims+2:4])
      ll = self$hlow
      ll = c(ll[Di], ll[self$dims+1], ll[2*self$dims+2:4])
      hp0 = topHP2
      hp0 = c(hp0[Di], hp0[self$dims+1], hp0[2*self$dims+2:4])
      dead = Optimizer_2(lhood,
                         grad_lhood,
                         N0    = 0,
                         x0    = hp0,
                         ran   = rbind(ll, uu),
                         maxevals = 100
      )
      
      cat("Fixed Wiggle Hyperparams optimized! ")
      return(topHP2)
      
    },
    HyperOptimise_fixedwiggle_finetune = function(with_prior=T){
      Di = 1:self$dims
    
      # First optimize a standard GP with i.i.d. noise
      dead  = self$Lhood_Prep()
   
      # Next make functions with one set of length scales for both GPs
      topLH  = -Inf
      topHP2  = 0
      lhood = function(hp){
        if(length(hp)!=self$dims+4)stop("finetune wiggle lhood wrong hp length")
        
        hp = c(hp[Di], hp[self$dims+1], hp[Di], hp[self$dims+2:4])
        tryCatch({
          O = self$Lhood_standard(hp) + with_prior*self$LHyperPrior(hp)
        }, error=function(e)-100000)
        if(O>topLH){
          topLH  <<- O
          topHP2 <<- hp
        }
        return(O)
      }
      grad_lhood = function(hp){
        if(length(hp)!=self$dims+4)stop("finetune wiggle grad lhood wrong hp length")
        
        hp = c(hp[Di], hp[self$dims+1], hp[Di], hp[self$dims+2:4])
        out = self$DLhood_standard(hp)+ with_prior*self$DLHyperPrior(hp)
        out[Di] = out[Di] + out[self$dims+1+Di]
        c(out[Di], out[self$dims+1], out[2*self$dims+2:4])
      }
      
      if(!is.null(self$HP)){
        hp0 = self$HP
        hp0 = c(hp0[Di], hp0[self$dims+1], hp0[2*self$dims+2:4])
        dead = lhood(hp0)
        topHP = self$HP
        bestV = sum(self$HP[self$dims*2 + 2:4])
      }else stop("cannot finetune fixed wiggle without last pars!")
      
      # Optimize the mixture of noise models
      topObj = -Inf
      topRho = 0
      Obj = function(rho){
        if(any(rho<0) | any(rho>1)) return(-1000000)
        x = c(topHP[1:(self$dims+1)],
              (1-rho[2])*rho[1]    *bestV,
              (1-rho[2])*(1-rho[1])*bestV,
              rho[2]*bestV)
        
        O = lhood(x)
        
        if(O>topObj){
          topObj <<- O
          topRho <<- rho
        }
        
        return( O )
      }
      dead = apply(LHCdim(40, 2), 1, Obj)
      dead = optim(f=Obj, par=topRho, method="Nelder-Mead", control = list(maxit=25, fnscale=-1))
      cat("Noise model mixture finetuned,...")

      # now optimize everything!
      uu = self$htop
      uu = c(uu[Di], uu[self$dims+1], uu[2*self$dims+2:4])
      ll = self$hlow
      ll = c(ll[Di], ll[self$dims+1], ll[2*self$dims+2:4])
      hp0 = topHP2
      hp0 = c(hp0[Di], hp0[self$dims+1], hp0[2*self$dims+2:4])
      dead = Optimizer_2(lhood,
                         grad_lhood,
                         N0    = 0,
                         x0    = hp0,
                         ran   = rbind(ll, uu),
                         maxevals = 100
      )
      
      cat("Fixed Wiggle Hyperparams finetuned! ")
      return(topHP2)
      
      
    },
    
    Refresh = function(Hpars=NULL, Ns=5, learnHpars=0, with_prior=T){
      
      if(self$copula){
        self$yd    = qnorm( order(self$yd_o)/(length(self$yd_o)+1) )
      }else{
        self$yd    = self$yd_o
      }
      
      xd         = self$xd
      yd         = self$yd
      
      
      self$ymean = mean(self$yd)
      self$Seeds = outer(xd[,self$dims+1], xd[,self$dims+1], "==")
      self$ydm   = self$yd - self$ymean
      self$SQ    = lapply(1:self$dims, function(d)outer(xd[,d], xd[,d], "-")^2)
      self$Noise = 0.0000*diag(length(yd))
      
      # learnHpars dictates what type of GP model we learn
      # 0: keep previous hyperparameters
      # 1: Compound Shperic, full optimisation
      # 2: compound Spheric, just fine tune update hpars
      # 3: wiggles
      # 4: wiggles fine tune
      # 5: unique seeds i.i.d. noise GP full optim
      # 6: unique seeds i.i.d. noise GP finetune
      # 7: use given hyperparameter, "Hpars"
      
      if(length(unique(self$xd[,self$dims+1]))==length(self$yd)){
        cat("no CRN data, doing unique seeds,")
        if(learnHpars==1) learnHpars=5
        if(learnHpars==2) learnHpars=6
        if(learnHpars==3) learnHpars=5
        if(learnHpars==4) learnHpars=6
      }
      
      if(learnHpars==0){
        cat(" keep old hyperpars, ")
        if (is.null(self$HP))stop("No old Pars! pass learnHpars=7 to use given Hpars")
        
      }else if(learnHpars==1){
        cat(" optimizing Hyperparameters, compound shperic, hyperprior ", with_prior, ", ")
        self$HP = self$HyperOptimise_CompSphere(with_prior=with_prior)
            
      }else if(learnHpars==2){
        cat(" finetuning Hyperparameters, compound spheric, hyperprior ", with_prior, ", ")
        self$HP = self$HyperOptimise_CompSphere_finetune(with_prior=with_prior)
            
      }else if(learnHpars==3){
        cat("optimising Hpers, fixed wiggles, prior, ")
        self$HP = self$HyperOptimise_fixedwiggle()
        
      }else if(learnHpars==4){
        cat("finetuning Hpers, fixed wiggles, prior, ")
        self$HP = self$HyperOptimise_fixedwiggle_finetune()
        
      }else if(learnHpars==5){
        cat(" optimizing Hpars, Unique seeds, ")
        self$HP = self$Unique_seeds_Optimise(with_prior=with_prior)
        
      }else if(learnHpars==6){
        cat(" finetuning Hpars, Unique seeds, ")
        self$HP = self$Unique_seeds_Optimise_finetune(with_prior=with_prior)
        
      }else if(learnHpars==7){
        cat(" using given hyperparameters, ")
        if (is.null(Hpars) | (length(Hpars)!=2*self$dims+4) | any(Hpars<0)) stop("cannot set hpars to NULL")
        self$HP = Hpars
      }
      
      V2 = self$HP[length(self$HP)]
      
      self$iK  = rcppeigen_invert_matrix( self$kernel(self$xd,self$xd) )
      self$iKY = self$iK%*%self$ydm
      self$Hpar_History[[length(self$yd)]] = self$HP
      
    },
    
    Optimise   = function(){
      topHP = exp(Optimizer(self$Lhood,
                            self$DLhood,
                            rbind(self$llow, self$ltop),
                            x0=log(c(LStheta,LSe,Sc)),
                            Ns=0
                            )$xmax)
      self$HP = topHP
      
      self$iK  = rcppeigen_invert_matrix(self$kernel(self$xd,self$xd) + 0.0001*diag(length(self$yd)))
      self$iKY = self$iK%*%(self$y-self$ydm)
      topHP
    },
    CDoptimise = function(){
      
      Freezer = function(fdims,fvals){
        odims = (1:self$npars)[!1:self$npars%in%fdims]
        initH0 = function()runif(1,self$hlow[odims],self$htop[odims])
        OBJ1 = function(input){ 
          LHP=rep(0,self$npars)
          for(i in 1:self$npars){
            if(i %in%fdims)LHP[i]=fvals[i==fdims] else LHP[i] = input[odims==i]} 
          self$Lhood(LHP)
        }
        DOBJ1 = function(input){ 
          LHP=rep(0,self$npars)
          for(i in 1:self$npars){if(i %in%fdims)LHP[i]=fvals[i==fdims] else LHP[i] = input[odims==i]} 
          self$DLhood(LHP)
        }
        c(initH0,OBJ1,DOBJ1)
      }
      
    },
    
    kernel_R = function(x1,x2){
      x1 = Check_X(x1, self$dims, T, "GP kernel x1")
      x2 = Check_X(x2, self$dims, T, "GP kernel x2")
      
      # return(Rcpp_Kernel2(x1,x2,self$HP))
      
      # cat(" x1 x2 dims ", dim(x1), " ", dim(x2)," ")
      
      Lt = self$HP[1:self$dims]
      St = self$HP[self$dims+1]
      Le = self$HP[self$dims+1 + 1:self$dims]
      Se = self$HP[self$dims+1 + 1+self$dims]
      Sc = self$HP[self$dims+1 + 1+self$dims +1]
      
      SQ = lapply(1:self$dims,function(d)outer(x1[,d],x2[,d],"-")^2)
      
      K_theta = St * exp(-0.5*Reduce("+",mapply("*", SQ,1/Lt^2,SIMPLIFY = F)))
      
      K_eps = Se * exp(-0.5*Reduce("+",mapply("*", SQ,1/Le^2,SIMPLIFY = F))) + Sc
      
      seed_mask = outer(x1[,ncol(x1)],x2[,ncol(x1)],"==")*(x1[,ncol(x1)]>0)
      
      output = K_theta + K_eps*seed_mask
      
      # cat("Cov matyrix dim ", dim(output),"\n") 
      output
      
      # browser()
      
    },
    kernel = function(x1,x2){
      # x1 = Check_X(x1, self$dims, T, "GP kernel x1")
      # x2 = Check_X(x2, self$dims, T, "GP kernel x2")
      
      return(Rcpp_Kernel2(x1,x2,self$HP))
    },
    dkernel.dx1_R = function(x1, x2){
      
      x1 = Check_X(x1, self$dims, T, "GP dkernel/dx1")
      x2 = Check_X(x2, self$dims, T, "GP dkernel/dx2")
      
      # DK = Rcpp_dKernel_dx1(x1,x2,self$HP)
      # return(lapply(1:self$dims, function(d)DK[1:nrow(x1) + nrow(x1)*(d-1), ,drop=F]))
      
      Lt = self$HP[1:self$dims]
      St = self$HP[self$dims+1]
      Le = self$HP[self$dims+1 + 1:self$dims]
      Se = self$HP[self$dims+1 + 1+self$dims]
      Sc = self$HP[self$dims+1 + 1+self$dims +1]
      
      DX = lapply(1:self$dims,function(d)outer(x1[,d],x2[,d],"-"))
      SQ = lapply(DX, function(dx)dx^2)
      
      K_theta = St * exp(-0.5*Reduce("+",mapply("*", SQ,1/Lt^2,SIMPLIFY = F)))
      dK_theta = lapply(1:self$dims, function(d)DX[[d]]*(1/Lt[d]^2)*K_theta)
      
      K_eps = Se * exp(-0.5*Reduce("+",mapply("*", SQ,1/Le^2,SIMPLIFY = F))) #+ Sc
      dK_eps = lapply(1:self$dims, function(d)DX[[d]]*(1/Le[d]^2)*K_eps)
      
      seed_mask = outer(x1[,ncol(x1)],x2[,ncol(x1)],"==")*(x1[,ncol(x1)]>0)
      
      output = lapply(1:self$dims,function(d)-dK_theta[[d]] - dK_eps[[d]]*seed_mask)
      
      output
      
    },
    dkernel.dx1 = function(x1,x2){
      # x1 = Check_X(x1, self$dims, T, "GP dkernel/dx1")
      # x2 = Check_X(x2, self$dims, T, "GP dkernel/dx2")
      
      return(lapply(1:self$dims, function(d)DK[1:nrow(x1) + nrow(x1)*(d-1), ,drop=F]))
    },
    MU     = function(x){
      x = Check_X(x, self$dims, T, "MU")
      # x = matrix(x)
      self$ymean+self$kernel(x,self$xd)%*%self$iKY
      },
    DTheta = function(xi){
      xi = matrix(xi, 1)
      xi = Check_X(xi, self$dims, F, fname = "grad theta")
      dKx = self$dkernel.dx1(cbind(xi,0), self$xd)
      lapply(dKx, function(dKx_i)dKx_i%*%self$iKY)
    },
    COV    = function(x1,x2){
      self$kernel(x1,x2) -self$kernel(x1,self$xd)%*%self$iK%*%self$kernel(self$xd,x2) },
    dCOV.dx1 = function(x1, x2, iKx2=NULL){
      if(is.null(iKx2))iKx2 = self$iK%*%self$kernel(self$xd, x2)
      dx2 = self$dkernel.dx1(x1, x2)
      dxd = self$dkernel.dx1(x1, self$xd)
      lapply(1:self$dims, function(d)dx2[[d]] - dxd[[d]]%*%iKx2)
    },
    dCOVxx.dx = function(x){
      # x = Check_X(x, self$dims, T, "dCOV_xx0")
      x=matrix(x,1)
      # if(nrow(x)!=1)stop("only can do one point in dCOVxx!")
      
      dKx = self$dkernel.dx1(x, self$xd)
      Kx = self$kernel(self$xd, x)
      prod = lapply(dKx, function(dk)-dk%*%self$iK%*%Kx)
      lapply(prod, function(p) (p+t(p)))
    },
    dCOVxx0.dx = function(x){
      x = Check_X(x, self$dims, T, "dCOV_xx0")
      if(nrow(x)!=1)stop("only can do one point in dCOVxx0!")

      Kx0   = self$kernel(matrix(c(x[1:self$dims],0),1), self$xd)
      dKx0  = self$dkernel.dx1(matrix(c(x[1:self$dims],0),1), self$xd)
      
      Kx    = self$kernel(self$xd, x)
      dKx   = lapply(self$dkernel.dx1(x, self$xd), t)
      
      lapply(1:self$dims, function(d)-dKx0[[d]]%*%self$iK%*%Kx - Kx0%*%self$iK%*%dKx[[d]])
    },
    SIGT_xr = function(x1, xr, iKxr=NULL){
      if(is.null(iKxr))iKxr = self$iK%*%self$kernel(self$xd, xr)
      
      # x1 = Check_X(x1, self$dims, T, "SIGT x1")
      x1 = matrix(x1,1)
      if(length(x1)!=self$dims+1)stop("only can do one point in sigt!")
      
      iCxx  = 1/sqrt(abs(self$COV(x1, x1)[1]))
      
      C_xr  = self$kernel(x1, xr) - self$kernel(x1,self$xd)%*%iKxr
      
      C_x0  = self$COV( c(x1[1:self$dims],0), x1)[1]
      
      c(C_xr, C_x0) * iCxx
      
    },
    dSIGT.dx_xr = function(x1, xr, iKxr=NULL){

      if(is.null(iKxr))iKxr = self$iK%*%self$kernel(self$xd, xr)
      
      x1=matrix(x1,1)
      
      # x1 = Check_X(x1, self$dims, T, "dSIGT x1")
      if(length(x1)!=self$dims+1)stop("only can do one point in dsigt/dx!")

      Cxx   = self$COV(x1, x1)[1]
      iC_xx = 1/Cxx
      dC_xx = sapply( self$dCOVxx.dx(x1), I)
      
      C_xr  = as.numeric( self$kernel(x1, self$xd)%*%iKxr )
      dC_xr = sapply( self$dCOV.dx1(x1,xr,iKxr), I)
      
      C_x0  = self$COV( matrix(c(x1[1:self$dims],0),1), x1)[1]
      dC_x0 = sapply( self$dCOVxx0.dx(x1), I)
      
      C_x = c(C_xr, C_x0)
      dC_x = rbind(dC_xr, dC_x0)
      
      # browser()
      
      sqrt(iC_xx) * (dC_x - 0.5*iC_xx*outer(C_x, dC_xx))
    },
    
    RecX = function(Ns=NULL, maxevals=50,oldrecx=NULL,...){
      if(is.null(Ns)){
        Ns = prod(sapply(1:self$dims, function(d)1.5*(self$XRAN[2,d]-self$XRAN[1,d])/self$HP[d]))
        Ns = min(Ns, 4000)
      }
      cat("optim RecX, Ns ",Ns, ".....")
      x0 = as.numeric(self$xd[which.max(self$yd), 1:self$dims])
      x0 = rbind(x0, oldrecx)

      fun = function(xi)self$MU(matrix(c(xi,0),1))
      dfun = function(xi)unlist(self$DTheta(xi))
      O = Optimizer_2(fun, dfun, self$XRAN, x0=x0, N0=Ns,...)$xmax
      
      cat("done, ")
      
      O
    },
    
    Plotting   = function(){
      if(self$dims!=1) return(0)
      
      Ym = self$MU(cbind(1:100, 0))
      se = 2*sqrt( sapply(1:100, function(x)self$COV(c(x,0),c(x,0))))
      
      plot(c(1, 100), range(c(Ym+se, Ym-se)),col="white")
      lines(1:100, Ym, type="l")
      points(self$xd[,1], self$yd, pch=19)
      
      NN = length(self$yd)
      points(self$xd[NN,1], self$yd[NN], pch=19, cex=2)
      
      lines(1:100, Ym + se, col="blue")
      lines(1:100, Ym - se, col="blue")
    }
  ))

}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
## ACQUISITION FUNCTIONS


KGCBcpp   = function(a, b){

  if(all(abs(b)<0.000001))return(0)
  
  big = abs(b)>0.000001
  if(any(!big)){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
  }
  
  nn=length(a)
  if(nn==1){
    return(0)
  }
  
  O=order(b)
  b=b[O]
  a=a[O]
  
  while (any(b[-1] == b[-nn])){
    Ri = which(b[-1] == b[-nn])[1]
    if(a[Ri] > a[Ri+1]){
      a=a[-(Ri+1)]
      b=b[-(Ri+1)]
    }else if(a[Ri] <= a[Ri+1]) {
      a=a[-Ri]
      b=b[-Ri]
    }
    nn=length(a)
  }
  
  if(nn==1){
    return(0)
  }
  
  O = KGCB_cpp(c(0,a),c(0,b))
  
  Z=c(-Inf,O[-1,2],Inf)
  I=O[,1]
  
  pC=pnorm(Z)
  dC=dnorm(Z)
  
  sum(  a[I]*(pC[-1]-pC[-length(pC)]) - b[I]*(dC[-1]-dC[-length(pC)])  ) - max(a)
}

Build_ref_X = function(GP, rounding=F){
  
  xd   = GP$xd
  XRAN = GP$XRAN
  
 
  
  xd = xd[,1:ncol(XRAN)]
  xd = Check_X(xd, ncol(XRAN), F, "Xr builder1")
  
  Xr = rbind(LHCran(nrow(xd), XRAN),
             xd + rnorm(length(xd)))
  
  
  Xr = sapply(1:GP$dims, function(d){Xr[,d] = pmax(Xr[,d], XRAN[1,d]); pmin(Xr[,d], XRAN[2,d])})
  
  if(rounding){
    Xr = round(Xr)
    
    repeats = rep(F, nrow(Xr))
    for(i in 1:(nrow(Xr)-1)) for (j in (i+1):nrow(Xr))if(all(Xr[i,]==Xr[j,]))repeats[i] = T
    
    while(any(repeats)){
      Xr[repeats,] = round(LHCran(sum(repeats), XRAN))
      repeats = rep(F, nrow(Xr))
      for(i in 1:(nrow(Xr)-1)) for (j in (i+1):nrow(Xr))if(all(Xr[i,]==Xr[j,]))repeats[i] = T
    }
  }
  
  Xr = cbind(Xr, 0)
  
  Check_ref_X(Xr, ncol(XRAN), T, "Xr builder2")
}

Check_ref_X = function(Xr, dims, withseeds, fname){
  
  Xr = Check_X(Xr, dims, withseeds, fname)
  
  NN = nrow(Xr)
  
  if(nrow(Xr)>1){
    repeats = sapply(1:(NN-1), function(i)any(sapply((i+1):NN, function(j)all(Xr[i,]==Xr[j,]))))
  
    repeats = c(repeats, F)
  
    Xr = Xr[!repeats,,drop=F]
    
  }
  
  Xr = Check_X(Xr, dims, withseeds, paste(fname, "reps removed"))
  
  Xr
  
}

KGCBfilter = function(a,b){
  
  if(all(abs(b)<0.000001))return(0)

  big = abs(b)>0.000001
  if(any(!big)){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
  }
  
  
  nn=length(a)
  if(nn==1){
    return(0)
  }
  
  O=order(b)
  b=b[O]
  a=a[O]
  #sort (a,b) by ascending b
  
  
  while (any(b[-1] == b[-nn])){
    Ri = which(b[-1] == b[-nn])[1]
    if(a[Ri] > a[Ri+1]){
      a=a[-(Ri+1)]
      b=b[-(Ri+1)]
    }else if(a[Ri] <= a[Ri+1]) {
      a=a[-Ri]
      b=b[-Ri]
    }
    nn=length(a)
  }
  
  if(nn==1){
    return(0)
  }
  
  A     = 1
  Alast = 1
  C=-Inf
  # print(c(N0(),"alast",Alast,"nn", nn))
  tt =1
  while(Alast<nn){
    s1     = Alast 
    sI     = (s1+1):nn
    CI     = -(a[s1]-a[sI])/(b[s1]-b[sI])
    
    bestsI = which.min(CI)
    Alast  = sI[bestsI]
    A      = c(A,Alast)
    C      = c(C,CI[bestsI])
    # print(c(tt, N0(),"Alast", Alast, "nn", nn, "bestS", bestsI, "as1", a[s1], "asI", a[sI]))
    tt = tt + 1
  }
  C=c(C,Inf)
  
  pC=pnorm(C)
  dC=dnorm(C)
  
  pC_diff = (pC[-1]-pC[-length(pC)])
  
  dC_diff = (dC[-1]-dC[-length(dC)])
  
  grad_a = rep(0, length(a))
  grad_b = rep(0, length(a))
  
  grad_a[A] = pC_diff
  
  sum(  a[A]*pC_diff - b[A]*dC_diff ) - max(a)
  
}

KGCBfilter_grad = function(a,b){
  a_old = a
  
  grad_a = rep(0, length(a))
  grad_a[which.max(a)] = -1
  
  grad_b = rep(0, length(a))
  
  if(all(abs(b)<0.000001))return(0)
  
  big = abs(b)>0.000001
  if(any(!big)>0){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
    new_order = c(which(big), which(a_old==a[length(a)]))
  }else{
    new_order = 1:length(a)
  }
  
  nn=length(a)
  
  if(nn==1){
    return(0)
  }
  
  
  O=order(b)
  b=b[O]
  a=a[O]
  new_order = new_order[O]
  #sort (a,b) by ascending b
  
  while (any(b[-1] == b[-nn])){
    Ri = which(b[-1] == b[-nn])[1]
    if(a[Ri] > a[Ri+1]){
      a=a[-(Ri+1)]
      b=b[-(Ri+1)]
      new_order = new_order[-(Ri+1)]
    }else if(a[Ri] <= a[Ri+1]) {
      a=a[-Ri]
      b=b[-Ri]
      new_order=new_order[-Ri]
    }
    nn=length(a)
  }
  
  
  if(nn==1){
    return(0)
  }
  A     = 1
  Alast = 1
  C=-Inf
  # print(c(N0(),"alast",Alast,"nn", nn))
  tt =1
  while(Alast<nn){
    s1     = Alast 
    sI     = (s1+1):nn
    CI     = -(a[s1]-a[sI])/(b[s1]-b[sI])
    
    bestsI = which.min(CI)
    Alast  = sI[bestsI]
    A      = c(A, Alast)
    C      = c(C, CI[bestsI])
    # print(c(tt, N0(),"Alast", Alast, "nn", nn, "bestS", bestsI, "as1", a[s1], "asI", a[sI]))
    tt = tt + 1
  }
  C=c(C,Inf)
  
  pC=pnorm(C)
  dC=dnorm(C)
  
  a_diff = c(0, a[A]) - c(a[A], 0)
  b_diff = c(0, b[A]) - c(b[A], 0)
  
  # sum(  a[A]*pC_diff - b[A]*dC_diff ) - max(a)
  KG_result = sum(  a_diff*pC - b_diff*dC ) - max(a)
  
  # g_1 = pC[2] #+ 2*C[1]*dC[1
  # print(A)
  # g_2 = pC[3] - pC[2] -1 # + 2* (C[2]*dC[2] - C[1]*dC[1])
  
  grad_a[new_order[A]] = grad_a[new_order[A]] + pC[-1] - pC[-length(pC)]
  
  dC_diff = dC[-1] - dC[-length(dC)]
  C[1] = 0
  C[length(C)] = 0
  
  grad_b[new_order[A]] = -dC[-1] + dC[-length(dC)]
  # browser()
  
  return(list(KG=KG_result, dmu=grad_a, dsig=grad_b))
  
}

UniformDesign_X = function(N0, ran, TestFun, double=0.5, rounding=F, Ns=NULL){
  
  # Sample 'double' of the x points twice, so 1 / (1+dbl) of budget is unique x
  Nx = ceiling(N0/(1+double))
  
  XX = LHCran(Nx, ran)
  
  dims = ncol(ran)
  
  XX = XX[c(sample(1:Nx), sample(1:Nx,N0-Nx)), ,drop=F]
  
  if(rounding)XX = round(XX)
  
  # browser()
  # OO = order(XX[,1])
  # XX = XX[OO,,drop=F]
  
  if(is.null(Ns)) SS = 1:N0 else SS = rep(1:Ns, len=N0)
  
  XX = cbind(XX, SS)
  
  for(i in 1:(N0-1)){
    repped = sapply((i+1):N0, function(j)all(XX[i,]==XX[j,]))
    while(any(repped)){
      XX[i,dims+1] = sample(which(1:Ns!=XX[i,dims+1]), size=1)
      repped = sapply((i+1):N0, function(j)all(XX[i,]==XX[j,]))
    }
  }
  
  
  XX
}

UniformDesign = function(N0, ran, TestFun, double=0.5, rounding=F, Ns=NULL){
  
  # Sample 'double' of the x points twice, so 1 / (1+dbl) of budget is unique x
  Nx = floor(N0/(1+double))
  
  XX = LHCran(Nx, ran)
  
  dims = ncol(ran)
  
  XX = XX[c(1:Nx, sample(1:Nx,N0-Nx)), ,drop=F]
  
  # browser()
  OO = order(XX[,1])
  XX = XX[OO,,drop=F]
  
  if(is.null(Ns)) SS = 1:N0 else SS = rep(1:Ns, len=N0)
  
  XX = cbind(XX, SS)
  
  if(rounding)XX = round(XX)
  
  YY  = TestFun(XX)
  
  if(is.null(Ns)){
    list(x=XX[,1:dims, drop=F], y=YY)
  }else{
    list(x=XX, y=YY)
  }
  
}

CRNUniformDesign = function(Nx, Ns, XRAN, TestFun, rounding=F){
  
  XX = LHCran(Nx, XRAN)
  XX   = XX[rep(1:Nx, Ns),,drop=F]

  XXs  = cbind(XX, rep(1:Ns, each=Nx))
  if(rounding)XXs = round(XXs)
  
  YY   = TestFun(XXs)
  list(x=XXs, y=YY)
}

CRNUniformDesign_ATO = function(N0, Ns, XRAN, TestFun, rounding=F){
  
  XX = LHCran(N0, XRAN)
  
  XXs  = cbind(XX, rep(1:Ns, len=N0))
  
  if(rounding)XXs = round(XXs)
  YY   = TestFun(XXs)
  list(x=XXs, y=YY)
}

New_CRNUniformDesign = function(N0, XRAN, TestFun, split=0.5, rounding=F){
  
  N_s = 3
  
  N_x_d = round(N0*split/N_s)
  N_x = N0 - N_x_d*N_s
  
  # Get a set of X points
  XX = LHCran(N_x, XRAN)
  
  # cherry pick some for all seed sampling
  d_i   = sample(N_x, N_x_d)
  XX_d  = XX[d_i, ,drop=F]
  XX_d  = XX_d[rep(1:N_x_d, N_s), ,drop=F]
  XX_d  = cbind(XX_d, rep(1:N_s, each=N_x_d))
  
  # get the remaining X for individual seeds
  s_i   = (1:N_x)[! d_i %in% (1:N_x)]
  XX_s  = XX[s_i, ,drop=F]
  XX_s  = cbind(XX_s, rep(1:N_s, len=nrow(XX_s)))
  
  
  # append x values and 
  XX   = rbind(XX_d, XX_s)
  
  if(rounding) XX = round(XX)
  
  YY   = TestFun(XX)
  list(x=XX, y=YY)
}





Make_CRNKG_grad = function(CRNGP, Xr){
  
  # Check reference points and remove repeats
  Xr = Check_ref_X(Xr, CRNGP$dims, T, "CRNKG Xr")
  ref_seed = Xr[1,CRNGP$dims+1]
  
  # browser()
  iKr = CRNGP$kernel(Xr, CRNGP$xd)%*%CRNGP$iK
  Mr  = CRNGP$MU(Xr)
  
  D = 1:CRNGP$dims
  tXr = t(Xr[D,])
  
  V2 = 0 #CRNGP$HP[length(CRNGP$HP)]
  
  CRNKG = function(xs){
    
    # repeats = apply(xs[D]==tXr, 2, all)
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats] }
    
    xs = matrix(xs,1)
    # xs = Check_X(xs, CRNGP$dims, T, "CRNKG xs")
    
    # if x is in the set Xr
    # repeats = apply(Xr, 1, function(xi)all(xi[D]==xs[D]))
    
    
    SDx = sqrt(abs(CRNGP$COV(xs, xs)) + V2 )[1]
    
    # browser()
    
    P0 = CRNGP$kernel(Xr, xs)
    # Nr = P0 > 0.001
    # 
    SIGT = rep(0, length(P0))
    # 
    # SIGT[Nr] = ( P0[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, xs) ) / SDx
    SIGT = ( P0 - iKr%*%CRNGP$kernel(CRNGP$xd, xs) ) * (1/SDx)
    
    xs0  = matrix(c(xs[D], ref_seed), 1)
    
    SIGT = c(SIGT, CRNGP$COV(xs0, xs)/SDx)
    
    MM   = c(Mr, CRNGP$MU(xs0))
    
    O = KGCBfilter(MM, SIGT)
    
    O
  }
  
  dCRNKG = function(xs){
    # xs = Check_X(xs, CRNGP$dims, T, "CRNKG xs")
    xs = matrix(xs, 1)
    # if x is in the set Xr
    # repeats = apply(Xr, 1, function(xi)all(xi[D]==xs[D]))
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats] }
    
    SDx = sqrt(abs(CRNGP$COV(xs, xs)) + V2)[1]
    
    # browser()
    
    P0 = CRNGP$kernel(Xr, xs)
    Nr = P0 > -Inf
    # 
    SIGT = rep(0, length(P0))
    # 
    SIGT[Nr] = ( P0[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, xs) ) / SDx
    # SIGT = ( P0 - iKr%*%CRNGP$kernel(CRNGP$xd, xs) ) * (1/SDx)
    
    xs0  = matrix(c(xs[D], ref_seed), 1)
    
    SIGT = c(SIGT, CRNGP$COV(xs0, xs)/SDx)
    
    MM   = c(Mr, CRNGP$MU(xs0))
    
    O = KGCBfilter_grad(MM, SIGT)
    
    Nr = abs(O$dsig)>0
    Nrr = which(Nr[-length(Nr)])
    
    Nr = which(Nr)
    
    # browser()
    
    dKG_s = (O$dsig[Nr]%*%CRNGP$dSIGT.dx_xr(xs, Xr[Nrr,,drop=F], t(iKr[Nrr,,drop=F]))  )[1:CRNGP$dims]
    
    dKG_m = unlist( CRNGP$DTheta(xs[1:CRNGP$dims])) * O$dmu[length(O$dmu)]
    # dKG_s = (O$dsig%*%CRNGP$dSIGT.dx_xr(xs, Xr, t(iKr)))[1:CRNGP$dims]
    
    list(dKG=dKG_m + dKG_s, KG=O$KG)
  }
  
  return(list(KG=CRNKG, dKG=dCRNKG))
}

MCMC_CRNKG_grad = function(CRNGPs, Xr, check_Seeds=NULL,
                           N0=1000, Na=1, maxevals=50){
  
  dims = CRNGPs[[1]]$dims
  
  if(is.null(check_Seeds)) check_Seeds = max(CRNGPs[[1]]$xd[,dims+1])
  
  cat("optim CRNKG_grad, check seeds:", check_Seeds, "...")
  
  # Check reference points and remove repeats
  Xr = Check_ref_X(Xr, CRNGPs[[1]]$dims, T, "CRNKG Xr")
  
  topKG = -Inf
  topxs = 0
  
  if(length(CRNGPs)==1){
    KG_FUNS = Make_CRNKG_grad(CRNGPs[[1]], Xr)
    KG_i = KG_FUNS$KG
    dKG_i = KG_FUNS$dKG
    KG = function(xs){
      KG_v = KG_i(xs)
      if(KG_v > topKG){
        
        topKG <<- KG_v
        topxs <<- xs

        # cat("best x ", signif(topxs,3), "  value ", signif(topKG),"\n")
      }
      KG_v
    }
    
    dKG = function(xs){
      O = dKG_i(xs)
      if(O$KG > topKG){
        topKG <<- O$KG
        topxs <<- xs
      }
      return(O$dKG)
    }
    
  }else{
    KG_FUNS = lapply(CRNGPs, function(CRNGP)Make_CRNKG_grad(CRNGP, Xr) )
    
    KG = function(xs){
      KG_v = mean(sapply(KG_FUNS, function(KG_F)KG_F$KG(xs)))
      if(KG_v > topKG){
        topKG <<- KG_v
        topxs <<- xs
      }
      KG_v
    }
    
    dKG = function(xs){
      KG_v = sapply(KG_FUNS, function(KG_F)KG_F$dKG(xs))
      KG_v = apply(KG_v, 1, mean)
      KG_v
    }
    
  }
  
  
  # KG_i = sapply(check_seeds, function(s)sapply(1:99,function(xi)KG(c(xi,s))))
  # plot(c(0,100), range(KG_i), col="white")
  # BB = sapply(check_seeds, function(s)lines(1:99, KG_i[,s], col=s+1))
  
  for(s in check_Seeds){
    KGs = function(x)KG(c(x,s))
    dKGs = function(x)dKG(c(x,s))
    dead = Optimizer_2(KGs, dKGs, ran=CRNGPs[[1]]$XRAN,
                       N0=N0,
                       Na=Na,
                       maxevals=maxevals)
  }
  
  # cat("best KG ", topxs, "\t", topKG,"\t", max(KG_i), "\n")
  cat("done, ")
  topxs
}




Make_PWKG_grad = function(CRNGP, Xr, ratio=1){
  Xr = Check_ref_X(Xr, CRNGP$dims, T, "PWKG Xr")
  new_seed = max(CRNGP$xd[,CRNGP$dims+1]) + 1
  
  iKr = CRNGP$kernel(Xr, CRNGP$xd)%*%CRNGP$iK
  Mr  = CRNGP$MU(Xr)
  
  D = 1:CRNGP$dims
  
  tXr = t(Xr[D,])
  
  KGFUNS = Make_CRNKG_grad(CRNGP, Xr)
  
  topKG = -Inf
  topX  = 0
  V2 = CRNGP$HP[length(CRNGP$HP)]
  
  KG = function(x){
    xs = matrix(c(x, new_seed), 1)
    # xs = Check_X(xs, CRNGP$dims, T, "PWKG KG xs")
    O = KGFUNS$KG(xs)
    
    if(O>topKG){
      topKG <<- O
      topX <<-  xs
    }
    
    O
  }
  
  dKG = function(x){
    xs = matrix(c(x, new_seed), 1)
    # xs = Check_X(xs, CRNGP$dims, T, "PWKG KG xs")
    
    O = KGFUNS$dKG(xs)
    
    if(O$KG>topKG){
      topKG <<- O$KG
      topX <<-  xs
    }
    
    return(O$dKG)
  }
  
  PWKG = function(xs1xs2){
    
    x1 = matrix( c(xs1xs2[D], new_seed) , 1)
    x2 = matrix( c(xs1xs2[CRNGP$dims + D], new_seed), 1)
    
    # repeats = apply(x1[D]==tXr, 2, all) | apply(x2[D]==tXr, 2, all)
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats] }
    
    # x1 = Check_X(x1, CRNGP$dims, T, "PWKG PW x1")
    # x2 = Check_X(x2, CRNGP$dims, T, "PWKG PW x2")
    
    if (all(x1==x2)) return(0)
    

    SDx = sqrt(abs(abs(CRNGP$COV(x1, x1)) + abs(CRNGP$COV(x2, x2)) - 2*CRNGP$COV(x1, x2)) + 2*V2)[1]
    
    P1 = CRNGP$kernel(Xr, x1)
    P2 = CRNGP$kernel(Xr, x2)
    Nr = P1 > 0.001 | P2 > 0.001
    
    SIGT = rep(0, length(P1))
    SIGT[Nr] =            P1[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x1) 
    SIGT[Nr] = SIGT[Nr] - P2[Nr] + iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x2)
    
    xs01  = matrix(c(x1[D], 0), 1)
    xs02  = matrix(c(x2[D], 0), 1)
    
    SIGTx1 = CRNGP$COV(xs01, x1)[1] - CRNGP$COV(xs01, x2)[1]
    SIGTx2 = CRNGP$COV(xs02, x1)[1] - CRNGP$COV(xs02, x2)[1]
    
    SIGT = c(SIGT, SIGTx1, SIGTx2) *(1 / SDx)
    
    MM   = c(Mr, CRNGP$MU(xs01), CRNGP$MU(xs02))
    
    # browser()
    
    O = 0.5 * KGCBfilter(MM, SIGT) * ratio
    
    if(O>topKG){
      topKG <<- O
      topX <<-  rbind(x1, x2)
    }
    
    O
    
  }
  
  dPWKG = function(xs1xs2){
    
    x1 = matrix( c(xs1xs2[D], new_seed), 1)
    x2 = matrix( c(xs1xs2[CRNGP$dims + D], new_seed), 1)
    
    # x1 = Check_X(x1, CRNGP$dims, T, "PWKG PW x1")
    # x2 = Check_X(x2, CRNGP$dims, T, "PWKG PW x2")
    
    if (all(x1==x2)) return(rep(0, CRNGP$dims*2))
    
    # if(xs[1]==75)browser()
    # if x is in the set Xr
    # repeats = apply(Xr, 1, function(xi)all(xi[D]==x1[D])) | apply(Xr, 1, function(xi)all(xi[D]==x2[D]))
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats]}
    # browser()
    VARx = abs(abs(CRNGP$COV(x1, x1)) + abs(CRNGP$COV(x2, x2)) - 2*CRNGP$COV(x1, x2) + 2*V2)[1]
    SDx  = sqrt(VARx)
    
    P1 = CRNGP$kernel(Xr, x1)
    P2 = CRNGP$kernel(Xr, x2)
    Nr = P1 > 0.001 | P2 > 0.001
    
    COV1     = rep(0, length(P1))
    COV1[Nr] = P1[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x1)
    
    COV2     = rep(0, length(P1))
    COV2[Nr] = P2[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x2)
    
    xs01  = matrix( c(x1[D], 0), 1)
    xs02  = matrix( c(x2[D], 0), 1)
    
    COV1 = c( COV1, CRNGP$COV(x1, xs01)[1], CRNGP$COV(x1, xs02)[1])
    COV2 = c( COV2, CRNGP$COV(x2, xs01)[1], CRNGP$COV(x2, xs02)[1])
    
    COV_dif = COV1 - COV2
    
    SIGT = COV_dif * (1 / SDx)
    
    MM   = c(Mr, CRNGP$MU(xs01), CRNGP$MU(xs02))
    
    O = lapply(KGCBfilter_grad(MM, SIGT), function(f)ratio*0.5*f)
    
    if(O$KG>topKG){
      topKG <<- O$KG
      topX <<-  rbind(x1,x2)
    }
    
    # browser()
    Nr = Nr & abs(O$dsig[1:(length(O$dsig)-2)])>0
    
    # browser()
    
    ul = function(L)sapply(L,I)
    
    # First dKG/dmu * dmu/dx12
    dkg.dmu = O$dmu[length(O$dmu) + rep(c(-1, 0),each=CRNGP$dims)]
    dMM = ul( c(CRNGP$DTheta(x1[D]), CRNGP$DTheta(x2[D])) )
    
    # Second dKG/dsig * dsig/dx12
    # sigt(x1, x2; Xr,x1,x2) = COV(x1, (Xr,x10,x20)) - COV(x2, (Xr,x10,x20)) / sqrt(SDx)
    # A=numerator dA/dx1 then d/dx2
    XrN = Xr[Nr,,drop=F]
    iKN = t(iKr[Nr,,drop=F])
    
    dCOV1 = matrix(0, nrow(Xr), CRNGP$dims)
    if(sum(Nr)>0) dCOV1[Nr,] = ul(CRNGP$dCOV.dx1(x1, XrN, iKN))
    dCOV1 = rbind(dCOV1, 
                  ul(CRNGP$dCOVxx0.dx(x1)) - ul(CRNGP$dCOV.dx1(xs01, x2)), 
                  ul(CRNGP$dCOV.dx1(x1, xs02)))
    
    
    dCOV2 = matrix(0, nrow(Xr), CRNGP$dims)
    if(sum(Nr)>0) dCOV2[Nr,] = ul(CRNGP$dCOV.dx1(x2, XrN, iKN))
    dCOV2 = rbind(dCOV2, 
                  ul(CRNGP$dCOV.dx1(x2, xs01)), 
                  ul(CRNGP$dCOVxx0.dx(x2)) - ul(CRNGP$dCOV.dx1(xs02, x1)))
    
    dCOV = cbind(dCOV1, -dCOV2)
    
    # B=denominator^2 d/dx COV(x1,x1) + COV(x2,x2) - 2 COV(x1,x2)
    
    dSD1 = ul(CRNGP$dCOVxx.dx(x1)) - 2*ul(CRNGP$dCOV.dx1(x1,x2))
    dSD2 = ul(CRNGP$dCOVxx.dx(x2)) - 2*ul(CRNGP$dCOV.dx1(x2,x1))
    
    dSD = c(dSD1, dSD2)
    
    # dsig/dx12 = (dA - dB A /2B) / sqrt(B)
    
    dSIGT = (dCOV - outer(COV_dif, dSD)*(0.5/VARx)) *(1/SDx)
    
    dPWKG = dkg.dmu*dMM + as.numeric(O$dsig%*%dSIGT)
    
    dPWKG
    
  }
  
  best_x = function(){
    topX
  }
  
  return(list(KG=KG, dKG=dKG, PWKG=PWKG, dPWKG=dPWKG, best_x=best_x))
}

MCMC_PWKG_grad = function(CRNGPs, Xr, ratio=1, 
                          N0=1000, Na=1, maxevals=50, 
                          PN0=4000, PNa=4, Pmaxevals=100){
  
  stopifnot(
    typeof(CRNGPs)=="list",
    ncol(Xr)==CRNGPs[[1]]$dims+1
  )
  
  cat("optim PWKG_grad with pairs, ...")
  topKG = -Inf
  topxs = 0
  topPW = FALSE
  
  
  new_seed = max(CRNGPs[[1]]$xd[,CRNGPs[[1]]$dims+1]) + 1
  
  if(length(CRNGPs)==1){
    
    FUNS = Make_PWKG_grad(CRNGPs[[1]], Xr, ratio=ratio)
    
    KG = function(xs){
      O = FUNS$KG(xs)
      if (O>topKG){
        topKG <<- O
        topxs <<- c(xs, new_seed)
        topPW <<- F
        # cat("best single ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dKG = function(xs){
      FUNS$dKG(xs)
    }
    
    
    PWKG = function(xs1xs2){
      O = FUNS$PWKG(xs1xs2)
      if (O>topKG){
        topKG <<- O
        topxs <<- cbind(matrix(xs1xs2, 2, byrow = T), new_seed)
        topPW <<- T
        # cat("best pair ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dPWKG = function(xs1xs2){
      FUNS$dPWKG(xs1xs2)
    }
    
  }else{
    
    FUNS = lapply(CRNGPs, function(CRNGP)Make_PWKG_grad(CRNGP, Xr) )
  
    KG = function(xs){
      O = mean(sapply(FUNS, function(FF)FF$KG(xs)))
      if (O>topKG){
        topKG <<- O
        topxs <<- c(xs, new_seed)
        topPW <<- F
        # cat("best single ", topxs, " value ", topKG, "\n")
      }
      O
  }
    dKG = function(xs){
      dKG_v = sapply(FUNS, function(FF)FF$dKG(xs))
      O = apply(dKG_v, 1, mean)
      O
    }
  
    PWKG = function(xs1xs2){
      O = mean(sapply(FUNS, function(FF)FF$PWKG(xs1xs2)))
      if (O>topKG){
        topKG <<- O
        topxs <<- cbind(matrix(xs1xs2, 2, byrow = T), new_seed)
        topPW <<- T
        # cat("best pair ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dPWKG = function(xs){
      dKG_v = sapply(FUNS, function(FF)FF$dPWKG(xs))
      O = apply(dKG_v, 1, mean)
      O
    }
    
  }
  # par(mfrow=c(1,1))
  # XX = 1:100
  # XXXX = cbind(XX, rep(XX, each=100))
  # KG_eval = sapply(XX, function(xi)KG(c(xi, new_seed)))
  # PWKG_eval = matrix(apply(XXXX, 1, PWKG), 100, 100)
  # contour(XX, XX, PWKG_eval)
  # PWKG_eval = apply(PWKG_eval, 1, max)
  # # plot(range(XX), range(c(KG_eval, PWKG_eval)), col="white")
  # scale = 100 / max(c(KG_eval, PWKG_eval))
  # lines(XX, scale*KG_eval, col="red")
  # lines(XX, scale*PWKG_eval, col="blue", lwd=2)
  
  # browser()
  
  # A = EI_optimizer(KG, CRNGPs[[1]]$XRAN)

  # A = EI_optimizer(PWKG, cbind(CRNGPs[[1]]$XRAN, CRNGPs[[1]]$XRAN))
  
  A = Optimizer_2(KG, dKG, CRNGPs[[1]]$XRAN,
                  N0=N0,
                  Na=Na,
                  maxevals=maxevals)

  A = Optimizer_2(PWKG, dPWKG, cbind(CRNGPs[[1]]$XRAN, CRNGPs[[1]]$XRAN),
                  N0 = PN0,
                  Na = PNa,
                  maxevals = Pmaxevals)

  
  topxs = FUNS$best_x()
  
  if (length(topxs)>CRNGPs[[1]]$dims+1){
    cat("best is pair, ")
  }else{
    cat("best is singleton, ")
  }
  
  return(topxs)
  
}




cat(" Loaded CRNMay2018.R\n\n")


