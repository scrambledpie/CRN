#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
NumericMatrix KGCB_cpp(NumericVector a,NumericVector b){
  int n = a.size()-1;
  int s = 2;
  double z = 0.0;
  
  NumericVector I(n);
  NumericVector Z(n);
  bool loopagain = true;
  I[1]=1; I[2]=2;
  Z[1]=-10000000000000; Z[2]=(a[1]-a[2])/(b[2]-b[1]);
  
  for(int i=3; i<n+1; i++){
    do{
      z=(a[i]-a[I[s]])/(b[I[s]]-b[i]);
      if(z<Z[s]){
        loopagain = true;
        s=s-1;
      }else{
        loopagain = false;
        s=s+1;
        I[s]=i;
        Z[s]=z;
      }
    }while(loopagain);
    
  }
  
  NumericMatrix out(s,2);
  for(int i =0; i<s; i++){out(i,0)=I[i+1]; out(i,1)=Z[i+1];}
  return out;
}


// [[Rcpp::export]]
NumericMatrix Rcpp_Kernel_old(NumericMatrix x1,
                          NumericMatrix x2,
                          NumericVector Hpars
){
  NumericMatrix K(x1.nrow(), x2.nrow());
  int d = x1.ncol()-1;
  double sqt;
  double sqe;
  NumericVector iLt(d);
  NumericVector iLe(d);
  NumericVector sq(d);
  
  for (int k=0; k<d; k++){
    iLt(k) = -0.5/(Hpars(k)*Hpars(k));
    iLe(k) = -0.5/(Hpars(k+d+1)*Hpars(k+d+1));
  }
  
  
  for(int i=0; i<x2.nrow(); i++ ){
    for(int j=0; j<x1.nrow(); j++ ){
      
      for(int k=0; k<d; k++){sq(k) = (x1(j,k)-x2(i,k)) * (x1(j,k)-x2(i,k));}
      
      sqt = 0.0;
      for(int k=0; k<d; k++){ sqt += sq(k) * iLt(k); }
      K(j,i) = Hpars(d) * exp(sqt);

      if ( abs(x1(j,d)-x2(i,d))<0.1 ){
        sqe=0.0;
        for(int k=0; k<d; k++){ sqe += sq(k) * iLe(k);}
        K(j,i) += Hpars(2*d+2) + Hpars(2*d+1) * exp(sqe);
      }
      
      
    }  
  }
  
  return K;
}

// [[Rcpp::export]]
NumericMatrix Rcpp_Kernel(NumericMatrix x1,
                          NumericMatrix x2,
                          NumericVector Hpars
){
  NumericMatrix K(x1.nrow(), x2.nrow());
  int d = x1.ncol()-1;
  double sqt;
  double sqe;
  NumericVector iLt(d);
  NumericVector iLe(d);
  NumericVector sq(d);
  
  for (int k=0; k<d; k++){
    iLt(k) = -0.5/(Hpars(k)*Hpars(k));
    iLe(k) = -0.5/(Hpars(k+d+1)*Hpars(k+d+1));
  }
  
  
  for(int i=0; i<x2.nrow(); i++ ){
    for(int j=0; j<x1.nrow(); j++ ){
      
      for(int k=0; k<d; k++){sq(k) = (x1(j,k)-x2(i,k)) * (x1(j,k)-x2(i,k));}
      
      sqt = 0.0;
      for(int k=0; k<d; k++){ sqt += sq(k) * iLt(k); }
      K(j,i) = Hpars(d) * exp(sqt);
      
      if ( abs(x1(j,d)-x2(i,d))<0.1 ){
        sqe=0.0;
        for(int k=0; k<d; k++){ sqe += sq(k) * iLe(k);}
        K(j,i) += Hpars(2*d+2) + Hpars(2*d+1) * exp(sqe);
        
        if( sum(sq)<0.00000000000001){K(j,i) += Hpars(2*d+3);}
      }
      
      
    }  
  }
  
  return K;
}

// [[Rcpp::export]]
NumericMatrix Rcpp_Kernel2(NumericMatrix x1,
                          NumericMatrix x2,
                          NumericVector Hpars
){
  NumericMatrix K(x1.nrow(), x2.nrow());
  int d = x1.ncol()-1;
  double sqt;
  double sqe;
  NumericVector iLt(d);
  NumericVector iLe(d);
  NumericVector sq(d);
  
  for (int k=0; k<d; k++){
    iLt(k) = -0.5/(Hpars(k)*Hpars(k));
    iLe(k) = -0.5/(Hpars(k+d+1)*Hpars(k+d+1));
  }
  
  
  for(int i=0; i<x1.nrow(); i++ ){
    for(int j=0; j<x2.nrow(); j++ ){
      
      for(int k=0; k<d; k++){sq(k) = (x1(i,k)-x2(j,k)) * (x1(i,k)-x2(j,k));}
      
      sqt = 0.0;
      for(int k=0; k<d; k++){ sqt += sq(k) * iLt(k); }
      K(i,j) = Hpars(d) * exp(sqt);
      
      if ( abs(x1(i,d)-x2(j,d))<0.1 ){
        sqe=0.0;
        for(int k=0; k<d; k++){ sqe += sq(k) * iLe(k);}
        K(i,j) += Hpars(2*d+2) + Hpars(2*d+1) * exp(sqe);
      }
      
      
    }  
  }
  return K;
}

// [[Rcpp::export]]
NumericMatrix Rcpp_dKernel_dx1(NumericMatrix x1,
                               NumericMatrix x2,
                               NumericVector Hpars
){
  NumericMatrix K(x1.nrow(), x2.nrow());
  int d = x1.ncol()-1;
  double sqt;
  double sqe;
  double kt;
  double ke;
  NumericVector iLt(d);
  NumericVector iLe(d);
  NumericVector sq(d);
  NumericVector dx(d);
  NumericMatrix dK(d*x1.nrow(), x2.nrow());
  
  for (int k=0; k<d; k++){
    iLt(k) = 1/(Hpars(k)*Hpars(k));
    iLe(k) = 1/(Hpars(k+d+1)*Hpars(k+d+1));
  }
  
  for(int i=0; i<x1.nrow(); i++ ){
    for(int j=0; j<x2.nrow(); j++ ){
      
      sqt = 0.0;
      for(int k=0; k<d; k++){ 
        dx(k) = x1(i,k)-x2(j,k); 
        sq(k) = dx(k)*dx(k); 
        sqt += sq(k) * iLt(k);
      }
      
      kt = Hpars(d) * exp(-0.5*sqt);
      
      for(int k=0; k<d; k++){
        dK(k*x1.nrow()+i, j) = -dx(k)*iLt(k)*kt;
      }
      
      
      if ( abs(x1(i,d)-x2(j,d))<0.1 ){
        sqe=0.0;
        for(int k=0; k<d; k++) sqe += sq(k) * iLe(k);
        
        ke = Hpars(2*d+1) * exp(-0.5*sqe);
        
        for(int k=0; k<d; k++){
          dK(k*x1.nrow()+i, j) -= dx(k)*iLe(k)*ke;
        }
      }
      
    }  
  }
  return dK;
}


