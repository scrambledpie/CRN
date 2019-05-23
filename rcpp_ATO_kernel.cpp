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

// [[Rcpp::export]]
double Simulate_ATO(NumericVector x, 
                    NumericMatrix Arrival,
                    NumericVector ProdTimeKey,
                    NumericVector ProdTimeNonKey,
                    const NumericMatrix Items,
                    const NumericMatrix Products
){
  
  int nProds = 5;
  int nItems = 8;
  int nKey   = 6;
  // int nNKey  = nItems-nKey; 
  double Profit = 0.0;
  NumericVector Orders(nProds);
  IntegerVector A(nProds);
  NumericVector Inventory = x;
  NumericMatrix ItemAvail(nItems,100);
  double PrevTime = 0.0; //PrevMinTime
  double Time = 0.0; //MinTime
  double warmup = 20.0;
  double Tmax = 70.0;
  int Prod; //Minprod
  int ptime_i =0;
  int pntime_i =0;
  int p=1;
  int q=1;
  int r=0;
  bool keyavail;
  double usedI;
  double BaseClock;
  
  
  // Initialise orders
  for(int i=0; i<nProds; i++){
    Orders(i) = Arrival(i,0);
    A(i) = 1;
  }
  
  //Initialise ItemAvail
  for (int i = 0;  i < nItems*100; i++) { ItemAvail[i] = -1.0; }
  // for (int i = 0;  i < nItems; i++) { ItemAvail(i,0) = 0; }
  int tt = 0;
  
  
  while(min(Orders)<Tmax){
    tt +=1;
    // Rcout << "\n Orders " << Orders << "  ";
    Time = min(Orders);
    Prod = which_min(Orders);
    Orders(Prod) = Orders(Prod) + Arrival(Prod, A(Prod));
    A(Prod) +=1;
    
    if(Time>Tmax)break;
    // Rcout << std::endl<< " Time, Prod " << Time << "   " << Prod << "  Profit ";
    
    // Update inventory arrived since last time and update items requested but not arrived.
    for (q=0; q<nItems; q++){
      r=0;
      while(ItemAvail(q,r)>0.0 & ItemAvail(q,r)<Time){ r+=1;}
      Inventory(q) += r;
      p=0;
      while(ItemAvail(q,p)>0.0){ ItemAvail(q,p) = ItemAvail(q,p+r); p+=1; }
    }
    
    // print out itemavail
    // Rcout << "ItemAvail" << std::endl;
    // for(int i=0; i<nItems; i++){
    //   for(int j=0; j<5; j++){ Rcout << ItemAvail(i,j) << " ";}
    //   Rcout << std::endl;
    // }
    
    // Rcout << " Inventory " << Inventory << std::endl;
    // Check if all items for this ordered prod are available.
    keyavail = true;
    for(int i=0; i<nKey; i++){
      if (Products(Prod, i+1) > Inventory(i)){keyavail=false; }
    }
    
    // Rcout << " keyavail " << keyavail << "  ";
    // make the sale, update stock, request more stock, calculate profit
    if(keyavail){
      
      // Decrease key items inventory and place replenishment orders for the amount of key items used
      for(int i=0; i<nKey; i++){
        usedI = Products(Prod, i+1);
        
        // if(tt==maxits) Rcout <<"here boy " <<  usedI << " ";
        
        if(usedI>0){
          
          if(Time>warmup){Profit += Items(i,0)*usedI;}
          
          Inventory(i) -= usedI;
          
          int ia = 0;
          BaseClock = -10.0;
          while(ItemAvail(i,ia)>0.0){BaseClock=ItemAvail(i,ia); ia+=1;}
          
          if(BaseClock<Time)BaseClock = Time;
          
          for(int j=0; j<usedI; j++){
            ItemAvail(i,ia+j) = BaseClock + Items(i,2) + Items(i,3)*ProdTimeKey(ptime_i);
            ptime_i +=1;
            BaseClock = ItemAvail(i,ia+j);
          }
        }
      }
      
      // Decrease nonkey items inventory and place replenishment orders for the amount of key items used
      for(int i=nKey; i<nItems; i++){
        
        usedI = Products(Prod, i+1);
        
        if(usedI>0 & usedI <= Inventory(i)){
          
          if(Time>warmup) Profit += Items(i,0)*usedI;
          
          Inventory(i) -= usedI;
          
          int ia = 0;
          BaseClock = -10.0;
          while(ItemAvail(i,ia)>0.0){BaseClock = ItemAvail(i,ia); ia+=1;}
          
          if(BaseClock<Time)BaseClock=Time;
          
          for(int j=0; j<usedI; j++){
            ItemAvail(i,ia+j) = BaseClock + Items(i,2) + Items(i,3)*ProdTimeNonKey(pntime_i);
            pntime_i +=1;
            BaseClock = ItemAvail(i, ia+j);
          }
        }
      }
      
      // Rcout << " Inventory " << Inventory << "  ";
      
      // Add on profit of each product sold and subtract holding cost of inventory
      
    }// end if keyavail
    
    // Rcout << "profit is " << Profit << std::endl << std::endl<< std::endl;
    if(Time>warmup){
      for(int i=0; i<nItems; i++){  Profit -= (Inventory(i)*Items(i,1)) * (Time-PrevTime);  }
    }
    
    // Rcout << Profit << std::endl ;
    PrevTime = Time;
    // Time = 100;
  }
  
  // Rcout << " Time " << Time << std::endl;
  // Rcout << " Prod " << Prod << std::endl;
  // Rcout << " Orders " << Orders << std::endl;
  // Rcout << " p " << ptime_i << std::endl;
  // Rcout << " q " << pntime_i << std::endl;
  // Rcout << " profit " << Profit << std::endl;
  // Rcout << "keyavail " << keyavail << std::endl;
  // Rcout << "Inventory " << Inventory << std::endl;
  // Rcout << "ItemAvail" << std::endl;
  // for(int i=0; i<nItems; i++){
  //   for(int j=0; j<5; j++){ Rcout << ItemAvail(i,j) << " ";}
  //   Rcout << std::endl;
  // }
  // Rcout << std::endl << std::endl << std::endl;
  
  
  return Profit/(Tmax-warmup);
}

// [[Rcpp::export]]
int Simulate_ATO2(NumericVector x, 
                     NumericMatrix Arrival,
                     NumericVector ProdTimeKey,
                     NumericVector ProdTimeNonKey,
                     const NumericMatrix Items,
                     const NumericMatrix Products
){
  
  int nProds = 5;
  int nItems = 8;
  int nKey   = 6;
  // int nNKey  = nItems-nKey; 
  double Profit = 0.0;
  NumericVector Orders(nProds);
  IntegerVector A(nProds);
  NumericVector Inventory = x;
  NumericMatrix ItemAvail(nItems,100);
  double PrevTime = 0.0; //PrevMinTime
  double Time = 0.0; //MinTime
  double warmup = 20.0;
  double Tmax = 70.0;
  int Prod; //Minprod
  int ptime_i =0;
  int pntime_i =0;
  int p=1;
  int q=1;
  int r=0;
  bool keyavail;
  double usedI;
  double BaseClock;
  
  
  // Initialise orders
  for(int i=0; i<nProds; i++){
    Orders(i) = Arrival(i,0);
    A(i) = 1;
  }
  
  //Initialise ItemAvail
  for (int i = 0;  i < nItems*100; i++) { ItemAvail[i] = -1.0; }
  // for (int i = 0;  i < nItems; i++) { ItemAvail(i,0) = 0; }
  int tt = 0;
  
  
  while(min(Orders)<Tmax){
    tt +=1;
    // Rcout << "\n Orders " << Orders << "  ";
    Time = min(Orders);
    Prod = which_min(Orders);
    Orders(Prod) = Orders(Prod) + Arrival(Prod, A(Prod));
    A(Prod) +=1;
    
    if(Time>Tmax)break;
    // Rcout << std::endl<< " Time, Prod " << Time << "   " << Prod << "  Profit ";
    
    // Update inventory arrived since last time and update items requested but not arrived.
    for (q=0; q<nItems; q++){
      r=0;
      while(ItemAvail(q,r)>0.0 & ItemAvail(q,r)<Time){ r+=1;}
      Inventory(q) += r;
      p=0;
      while(ItemAvail(q,p)>0.0){ ItemAvail(q,p) = ItemAvail(q,p+r); p+=1; }
    }
    
    // print out itemavail
    // Rcout << "ItemAvail" << std::endl;
    // for(int i=0; i<nItems; i++){
    //   for(int j=0; j<5; j++){ Rcout << ItemAvail(i,j) << " ";}
    //   Rcout << std::endl;
    // }
    
    // Rcout << " Inventory " << Inventory << std::endl;
    // Check if all items for this ordered prod are available.
    keyavail = true;
    for(int i=0; i<nKey; i++){
      if (Products(Prod, i+1) > Inventory(i)){keyavail=false; }
    }
    
    // Rcout << " keyavail " << keyavail << "  ";
    // make the sale, update stock, request more stock, calculate profit
    if(keyavail){
      
      // Decrease key items inventory and place replenishment orders for the amount of key items used
      for(int i=0; i<nKey; i++){
        usedI = Products(Prod, i+1);
        
        // if(tt==maxits) Rcout <<"here boy " <<  usedI << " ";
        
        if(usedI>0){
          
          if(Time>warmup){Profit += Items(i,0)*usedI;}
          
          Inventory(i) -= usedI;
          
          int ia = 0;
          BaseClock = -10.0;
          while(ItemAvail(i,ia)>0.0){BaseClock=ItemAvail(i,ia); ia+=1;}
          
          if(BaseClock<Time)BaseClock = Time;
          
          for(int j=0; j<usedI; j++){
            ItemAvail(i,ia+j) = BaseClock + Items(i,2) + Items(i,3)*ProdTimeKey(ptime_i);
            ptime_i +=1;
            BaseClock = ItemAvail(i,ia+j);
          }
        }
      }
      
      // Decrease nonkey items inventory and place replenishment orders for the amount of key items used
      for(int i=nKey; i<nItems; i++){
        
        usedI = Products(Prod, i+1);
        
        if(usedI>0 & usedI <= Inventory(i)){
          
          if(Time>warmup) Profit += Items(i,0)*usedI;
          
          Inventory(i) -= usedI;
          
          int ia = 0;
          BaseClock = -10.0;
          while(ItemAvail(i,ia)>0.0){BaseClock = ItemAvail(i,ia); ia+=1;}
          
          if(BaseClock<Time)BaseClock=Time;
          
          for(int j=0; j<usedI; j++){
            ItemAvail(i,ia+j) = BaseClock + Items(i,2) + Items(i,3)*ProdTimeNonKey(pntime_i);
            pntime_i +=1;
            BaseClock = ItemAvail(i, ia+j);
          }
        }
      }
      
      // Rcout << " Inventory " << Inventory << "  ";
      
      // Add on profit of each product sold and subtract holding cost of inventory
      
    }// end if keyavail
    
    // Rcout << "profit is " << Profit << std::endl << std::endl<< std::endl;
    if(Time>warmup){
      for(int i=0; i<nItems; i++){  Profit -= (Inventory(i)*Items(i,1)) * (Time-PrevTime);  }
    }
    
    // Rcout << Profit << std::endl ;
    PrevTime = Time;
    // Time = 100;
  }
  
  // Rcout << " Time " << Time << std::endl;
  // Rcout << " Prod " << Prod << std::endl;
  // Rcout << " Orders " << Orders << std::endl;
  // Rcout << " p " << ptime_i << std::endl;
  // Rcout << " q " << pntime_i << std::endl;
  // Rcout << " profit " << Profit << std::endl;
  // Rcout << "keyavail " << keyavail << std::endl;
  // Rcout << "Inventory " << Inventory << std::endl;
  // Rcout << "ItemAvail" << std::endl;
  // for(int i=0; i<nItems; i++){
  //   for(int j=0; j<5; j++){ Rcout << ItemAvail(i,j) << " ";}
  //   Rcout << std::endl;
  // }
  // Rcout << std::endl << std::endl << std::endl;
  
  
  return ptime_i;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/* R
set.seed(1)
print(KGCB_cpp(rnorm(10),rnorm(10)))

*/

