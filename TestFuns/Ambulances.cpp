#include <Rcpp.h>
using namespace Rcpp;


double abs2(double x){
  return std::abs(x);
}


double Dist2Call(double xcall,
                 double ycall,
                 double xb,
                 double yb,
                 double xlc,
                 double ylc,
                 double driven){
  
  // If the ambulace is still attending a call
  if (driven<0){Rcout << "negative driven!!"<<std::endl; return(100000);}
  
  // Otherwise compute it's distance to the call
  double xcur, ycur, driven_x;
  
  // Is it already at base?
  if(driven > abs2(xlc-xb) + abs2(ylc-yb) ){
    xcur = xb;
    ycur = yb;
    
  // Where is it on the x dirn way to base?
  }else if (driven > abs2(ylc-yb)){
    driven_x = driven - abs2(ylc-yb);
    if (xlc>xb){
      xcur = xlc - driven_x;
    }else{
      xcur = xlc + driven_x;
    }
    ycur = yb;
    
  // Where is it on the way to base? Vertical travel first
  }else{
    if (ylc>yb){
      ycur = ylc - driven;
    }else{
      ycur = ylc + driven;
    }
    xcur = xlc;
  }
  
  return(abs2(ycur-ycall) + abs2(xcur-xcall));
  
}


// [[Rcpp::export]]
double Ambulances_Square(NumericMatrix bases, 
                             NumericVector CallTimes,
                             NumericMatrix CallLocs,
                             NumericVector ServTimes){
  

  double velfk = 60; double vf = velfk/30; double ivf = 1/vf;
  double velsk = 40; double vs = velsk/30;
  int nA = bases.nrow();
  int NumCalls = CallLocs.nrow();
  int NumLate = 0;

  double xcall, ycall, xb, yb, xlc, ylc;
  double Ddriven, Tgap, minD2call, D2call, depart;
  int closestA;
  
  NumericVector ExitTimes(NumCalls);
  NumericVector AmbArrTimes(NumCalls);
  NumericVector closestAhist(NumCalls);
  NumericVector minDcall(NumCalls);
  
  
  // Stores all info about ambulances
  NumericMatrix Ambulances(nA, 5);
  for (int j=0; j<nA; j++){
    Ambulances(j,1) = bases(j,0);
    Ambulances(j,2) = bases(j,1);
    Ambulances(j,3) = bases(j,0);
    Ambulances(j,4) = bases(j,1);
  }

  for(int t=0; t<NumCalls; t++){
    
    // Rcout << t << std::endl;

    // Oh no, someone is hurt! Where are they!?
    xcall = CallLocs(t,0);
    ycall = CallLocs(t,1);

    closestA = -1;
    minD2call = 1000;

    // Check all the ambulances! Are they free??
    for (int j=0; j<nA; j++){
      

      Tgap = CallTimes(t) - Ambulances(j,0);
      
      // Is this ambulance available?
      if(Tgap>0){
        Ddriven = Tgap*vs;
  
        xlc = Ambulances(j, 1);
        ylc = Ambulances(j, 2);
        xb  = Ambulances(j, 3);
        yb  = Ambulances(j, 4);
  
        D2call = Dist2Call(xcall, ycall, xb, yb, xlc, ylc, Ddriven);
        
        if(D2call<minD2call){
          minD2call = D2call;
          closestA = j;
        }
      }
      

    }

    // finished checking ambulances!

    if (closestA >-1){
      // If one was available, it departs ASAP!
      depart = CallTimes(t);

    }else{
      // No ambulance available, so let's get the quickest one!
      NumLate +=1;
      
      depart = min(Ambulances(_,0));
      closestA = which_min(Ambulances(_,0));
      
      xlc = Ambulances(closestA, 1);
      ylc = Ambulances(closestA, 2);
      
      minD2call = abs2(xlc-xcall) + abs2(ylc-ycall);
      
      // Rcout << "No free ambulances! Nearest one is free" << depart << std::endl;
    }

    AmbArrTimes(t)          = depart +  minD2call*ivf;
    ExitTimes(t)            = AmbArrTimes(t) + ServTimes(t);
    Ambulances(closestA, 0) = ExitTimes(t);
    Ambulances(closestA, 1) = xcall;
    Ambulances(closestA, 2) = ycall;

    closestAhist(t) = closestA;
    minDcall(t) = minD2call;
  }

  // return Rcpp::List::create(
  //   Rcpp::Named("AmbArr") = AmbArrTimes,
  //   Rcpp::Named("ExitTimes") = ExitTimes,
  //   Rcpp::Named("closestA") = closestAhist,
  //   Rcpp::Named("minDcall") = minDcall,
  //   Rcpp::Named("NumLate") = NumLate,
  //   Rcpp::Named("A") = Ambulances
  // );
  return(mean(AmbArrTimes-CallTimes));
}

// [[Rcpp::export]]
NumericVector Ambulances_Square_timed(NumericMatrix bases, 
                               NumericVector CallTimes,
                               NumericMatrix CallLocs,
                               NumericVector ServTimes,
                               double Simtime){
  
  
  double velfk = 60; double vf = velfk/30; double ivf = 1/vf;
  double velsk = 40; double vs = velsk/30;
  int nA = bases.nrow();
  int NumCalls = CallLocs.nrow();
  int NumLate = 0;
  
  double xcall, ycall, xb, yb, xlc, ylc;
  double Ddriven, Tgap, minD2call, D2call, depart;
  int closestA;
  
  NumericVector ExitTimes(NumCalls);
  NumericVector AmbArrTimes(NumCalls);
  NumericVector closestAhist(NumCalls);
  NumericVector minDcall(NumCalls);
  
  
  // Stores all info about ambulances
  NumericMatrix Ambulances(nA, 5);
  for (int j=0; j<nA; j++){
    Ambulances(j,1) = bases(j,0);
    Ambulances(j,2) = bases(j,1);
    Ambulances(j,3) = bases(j,0);
    Ambulances(j,4) = bases(j,1);
  }
  
  for(int t=0; t<NumCalls; t++){
    
    
    if (CallTimes(t)>Simtime){
      // End the simulation if simtime has elapsed
      
      // Rcout << t << std::endl;
      // Rcout << CallTimes(t-1) << std::endl;
      // Rcout << CallTimes(t) << std::endl;

      // double traveltime = 0.0;
      // for (int ti=0; ti<t; ti++){traveltime += AmbArrTimes(ti) - CallTimes(ti);};
      // return(traveltime);
      
      NumericVector JourneyTimes(t);
      for (int ti=0; ti<t; ti++){JourneyTimes(ti) = AmbArrTimes(ti) - CallTimes(ti);};
      
      // std::sort(JourneyTimes.begin(), JourneyTimes.end());
        
      return(AmbArrTimes - CallTimes);
      // return(mean(AmbArrTimes[CallTimes<Simtime]-CallTimes[CallTimes<Simtime]));
      
    }
    
    
    // Oh no, someone is hurt! Where are they!?
    xcall = CallLocs(t,0);
    ycall = CallLocs(t,1);
    
    closestA = -1;
    minD2call = 1000;
    
    // Check all the ambulances! Are they free??
    for (int j=0; j<nA; j++){
      
      
      Tgap = CallTimes(t) - Ambulances(j,0);
      
      // Is this ambulance available?
      if(Tgap>0){
        Ddriven = Tgap*vs;
        
        xlc = Ambulances(j, 1);
        ylc = Ambulances(j, 2);
        xb  = Ambulances(j, 3);
        yb  = Ambulances(j, 4);
        
        D2call = Dist2Call(xcall, ycall, xb, yb, xlc, ylc, Ddriven);
        
        if(D2call<minD2call){
          minD2call = D2call;
          closestA = j;
        }
      }
      
      
    }
    
    // finished checking ambulances!
    
    if (closestA >-1){
      // If one was available, it departs ASAP!
      depart = CallTimes(t);
      
    }else{
      // No ambulance available, so let's get the quickest one!
      NumLate +=1;
      
      depart = min(Ambulances(_,0));
      closestA = which_min(Ambulances(_,0));
      
      xlc = Ambulances(closestA, 1);
      ylc = Ambulances(closestA, 2);
      
      minD2call = abs2(xlc-xcall) + abs2(ylc-ycall);
      
      // Rcout << "No free ambulances! Nearest one is free" << depart << std::endl;
    }
    
    AmbArrTimes(t)          = depart +  minD2call*ivf;
    ExitTimes(t)            = AmbArrTimes(t) + ServTimes(t);
    Ambulances(closestA, 0) = ExitTimes(t);
    Ambulances(closestA, 1) = xcall;
    Ambulances(closestA, 2) = ycall;
    
    closestAhist(t) = closestA;
    minDcall(t) = minD2call;
  }
  
  // return Rcpp::List::create(
  //   Rcpp::Named("AmbArr") = AmbArrTimes,
  //   Rcpp::Named("ExitTimes") = ExitTimes,
  //   Rcpp::Named("closestA") = closestAhist,
  //   Rcpp::Named("minDcall") = minDcall,
  //   Rcpp::Named("NumLate") = NumLate,
  //   Rcpp::Named("A") = Ambulances
  // );
  return(AmbArrTimes-CallTimes);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/** R
set.seed(1)
runlength = 100
NumCalls = runlength*30

Arrivals = rexp(NumCalls, rate = 1/60)
Arrivals = cumsum(Arrivals)

CallLocs = matrix(0, NumCalls, 2)
i=1
while (CallLocs[NumCalls, 1]==0){
  u = runif(3)
  if(1.6*u[3] <= 1.6 - abs(u[1]-0.8) - abs(u[2]-0.8)){
    CallLocs[i,] = u[1:2]
    i = i+1
  }
}


Stimes   = rchisq(NumCalls, 3)

Xi = c(0.25, 0.25, 0.75, 0.75)
Yi = c(0.25, 0.75, 0.25, 0.75)

X = cbind(Xi, Yi)

X = matrix(runif(8), 4,2)

Arrt = Ambulances_Square(X, Arrivals, CallLocs, Stimes)
*/
