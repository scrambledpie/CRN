#include <Rcpp.h>
using namespace Rcpp;

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
  int kk = 0;
  
  
  while(min(Orders)<Tmax){
    // Rcout << "\n Orders " << tt  << " "<< Orders << "  ";
    kk +=1;
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
      
      // Rcout << "\n Profit " << Profit  << " "<< kk << "  ";
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
      
      // Rcout << "\n Profit " << Profit  << " "<< kk << "  ";
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



