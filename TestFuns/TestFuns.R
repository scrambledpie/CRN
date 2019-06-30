library(Rcpp)

if (F){
  Check_X = function(x, dims, withseeds, fname){
    # Args
    #   x: input vec or matrix
    #   dims: number of continuous dims of x
    #   withseeds: should the x include an integer seed aslast element
    #   fname: printed with the error message
    #
    # Returns
    #   x: correctly formatted as a matrix
    
    dimsN = dims + withseeds
    
    if (dims==1&is.null(dim(x))&!withseeds){ cat(fname, "converting x to matrix"); x = matrix(x, ncol=1)}
    
    if (length(x)==dimsN & length(dim(x))==0){ cat(fname, "converting x to matrix"); x = matrix(x, ncol=dimsN, nrow=1)}
    
    if (length(x)%%dimsN!=0) stop(fname, " wrong length input, got ", length(x))
    
    if (length(dim(x))!=2) stop(fname, "wrong input dim, got ", dim(x))
    
    if (ncol(x)!=dimsN) stop(fname, "wrong # input columns, got ", ncol(x))
    
    if (withseeds & !all(x[,dimsN]%%1==0)) stop(fname, "non-integer seeds not allowed")
    
    if (withseeds & !all(x[,dimsN]>=0)) stop(fname, "negative seeds not allowed")
    
    x
  }
}else{
  Check_X = function(x, dims, withseeds, fname) x
}

LHCdim = function(N0, dims){
  # Args
  #   N0: the number of points to generate
  #   dims: the box for scaling+centring samples
  #   
  # Returns
  #   a randomly generated LHC of N0 points in [0,1]^dims
  
  sapply(1:dims, function(d){  (sample(1:N0) - runif(N0))  })*(1/N0)
}

LHCran = function(N0, ran){
  # Args
  #   N0: the number of points to generate
  #   ran: the box for scaling+centring samples
  #   
  # Returns
  #   a randomly generated LHC of N0 points in box ran dimensions
  
  if(length(ran)==2)ran=matrix(ran,ncol=1)
  
  if (any(ran[1,]>=ran[2,])|nrow(ran)!=2) stop("LHC bounds are not valid: ", paste(bottoms, tops, collapse=", "))
  
  dims = ncol(ran)
  LHC1 = sapply(1:dims,function(d){  (sample(1:N0) - runif(N0))  })*(1/N0)
  
  BB = matrix(ran[1,], N0, length(ran[1,]), byrow=T)
  CC = matrix(ran[2,]-ran[1,], N0, length(ran[1,]), byrow=T)
  
  BB + LHC1*CC
}

########################################################################################################
# ATO

old_RunATO = function(BaseStockLevel=rep(1,8), seed=1, runlength=5){

  # Args:
  #   BaseStockLEvel: target stock quantities for each product
  #   runlength: integer, repetitions of the simulation
  #   seed: integer, to force deterministic simulation
  #
  # Returns:
  #   fnAvg: Average profit over runlength reps of simulation
  
  NumComponentType=8;               #% number of component types
  
  ProTimeMean = c(0.15, 0.4, 0.25, 0.15, 0.25, 0.08, 0.13, 0.4)
  ProTimeStd = 0.15*ProTimeMean
  Profit = 1:8
  HoldingCost = matrix(2, 1, NumComponentType)
  # % holding cost of each component in
  # % inventory
  # 
  # % Parameters of Customers
  
  ArrivalRate = 12             #% assume Possion arrival
  # %NumCustomerType=5;             % number of customer types
  # %CustomerProb=[0.3,0.25,0.2,0.15,0.1];   % probability of each customer
  KeyComponent=c(1,0,0,1,0,1,0,0,
                 1,0,0,0,1,1,0,0,
                 0,1,0,1,0,1,0,0,
                 0,0,1,1,0,1,0,0,
                 0,0,1,0,1,1,0,0)
  KeyComponent = matrix(KeyComponent, nrow=5, byrow=TRUE)
  
  NonKeyComponent=c(0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,1,
                    0,0,0,0,0,0,1,0)
  NonKeyComponent = matrix(NonKeyComponent, nrow=5, byrow=TRUE)
  
  WarmUp = 20
  TotalTime = 70
  
  #% Simulation
  fnSum = 0.0
  set.seed(seed)
  
  for (k in 1:runlength){
    
    # runif_preset = runif()
    
    EventTime = matrix(1e5, 1, 1+NumComponentType)
    EventTime[1] = -log(runif(1))/ArrivalRate
    TotalProfit = 0
    TotalCost = 0
    Inventory = BaseStockLevel
    
    Clock = 0
    
    while (Clock<TotalTime){
      OldInventory = Inventory
      OldClock = Clock
      Clock = min(EventTime)
      event = which.min(EventTime)
      if (event==1){ # a customer has come!
        
        temp=runif(1)
        if (temp<0.3){
          CustomerType=1
        }else if (temp<0.55){
          CustomerType=2
        }else if (temp<0.75){
          CustomerType=3
        }else if (temp<0.9){
          CustomerType=4
        }else{
          CustomerType=5
        }
        
        
        
        KeyOrder = KeyComponent[CustomerType,]
        NonKeyOrder = NonKeyComponent[CustomerType,]
        Sell=1
        for (i in 1:NumComponentType){
          if (Inventory[i] - KeyOrder[i] < 0){
            Sell=0 
          }
          if (Inventory[i] - NonKeyOrder[i] < 0) {
            NonKeyOrder[i] = Inventory[i]
          }
        }
        if (Sell==1){
          Inventory = Inventory - KeyOrder - NonKeyOrder
          if (Clock>WarmUp){
            TotalProfit = TotalProfit + sum(Profit*(KeyOrder + NonKeyOrder))
          }
        }
        
        EventTime[1]=Clock-log(runif(1))/ArrivalRate
        if (Sell==1){
          for (i in 1:NumComponentType){
            if (Inventory[i]<BaseStockLevel[i]  && EventTime[i+1]>1e4){
              EventTime[i+1] = Clock + max(0,ProTimeMean[i]+rnorm(1)*ProTimeStd[i])
            }
          }
        }
      }else{ # stock has arrived
        ComponentType = event-1;
        Inventory[ComponentType] = Inventory[ComponentType]+1
        if (Inventory[ComponentType]>=BaseStockLevel[ComponentType]){
          EventTime[event]=1e5
          if (Clock>WarmUp){
            TotalCost=TotalCost+(Clock-OldClock)* sum(OldInventory*HoldingCost)
          }
        }
      }
    }
    fn = (TotalProfit-TotalCost)/(TotalTime-WarmUp);
    fnSum = fnSum + fn;
  }
  
  fnAvg = fnSum / runlength
  
  return(fnAvg)
}

RunATO = function(BaseStockLevel=rep(1,8), seed=1, runlength=5){

  # Args:
  #   BaseStockLEvel: target stock quantities for each product
  #   runlength: integer, repetitions of the simulation
  #   seed: integer, to force deterministic simulation
  #
  # Returns:
  #   fnAvg: Average profit over runlength reps of simulation
  
  NumComponentType=8;               #% number of component types
  
  ProTimeMean = c(0.15, 0.4, 0.25, 0.15, 0.25, 0.08, 0.13, 0.4)
  ProTimeStd = 0.15*ProTimeMean
  Profit = 1:8
  HoldingCost = matrix(2, 1, NumComponentType)
  # % holding cost of each component in
  # % inventory
  # 
  # % Parameters of Customers
  
  ArrivalRate = 12             #% assume Possion arrival
  iArrivalRate = 1/ArrivalRate
  # %NumCustomerType=5;             % number of customer types
  # %CustomerProb=[0.3,0.25,0.2,0.15,0.1];   % probability of each customer
  KeyComponent=c(1,0,0,1,0,1,0,0,
                 1,0,0,0,1,1,0,0,
                 0,1,0,1,0,1,0,0,
                 0,0,1,1,0,1,0,0,
                 0,0,1,0,1,1,0,0)
  KeyComponent = matrix(KeyComponent, nrow=5, byrow=TRUE)
  
  NonKeyComponent=c(0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,1,
                    0,0,0,0,0,0,1,0)
  NonKeyComponent = matrix(NonKeyComponent, nrow=5, byrow=TRUE)
  
  WarmUp = 20
  TotalTime = 70
  
  #% Simulation
  fnSum = 0.0
  set.seed(seed)
  Customer_series = sample(size=TotalTime*ArrivalRate*runlength*2, 1:5, prob=c(0.3,0.25,0.2,0.15,0.1), replace=T)
  t = 1
  
  for (k in 1:runlength){
    
    # runif_preset = runif()
    
    EventTime = matrix(1e5, 1, 1+NumComponentType)
    EventTime[1] = -log(runif(1))*iArrivalRate
    TotalProfit = 0
    TotalCost = 0
    Inventory = BaseStockLevel
    
    Clock = 0
    
    while (Clock<TotalTime){
      OldInventory = Inventory
      OldClock = Clock
      Clock = min(EventTime)
      event = which.min(EventTime)
      if (event==1){ # a customer has come!
        # Reset clock
        EventTime[1]=Clock-log(runif(1))*iArrivalRate
        
        # Get new customer
        CustomerType = Customer_series[t]
        t=t+1
        
        KeyOrder = KeyComponent[CustomerType,]
        NonKeyOrder = NonKeyComponent[CustomerType,]
        
        Sell = all(Inventory>=KeyOrder)
        NonKeyOrder = Inventory * (Inventory<=NonKeyOrder) + NonKeyOrder*(Inventory>NonKeyOrder)
        
        if (Sell){
          Inventory = Inventory - KeyOrder - NonKeyOrder
          if (Clock>WarmUp){ TotalProfit = TotalProfit + sum(Profit*(KeyOrder + NonKeyOrder)) }
          
          for (i in 1:NumComponentType){
            if (Inventory[i]<BaseStockLevel[i]  && EventTime[i+1]>1e4){
              EventTime[i+1] = Clock + max(0,ProTimeMean[i]+rnorm(1)*ProTimeStd[i])
            }
          }
          # reset = ((Inventory<BaseStockLevel)  & (EventTime[-1]>1e4))
          # new_times = rnorm(sum(reset), ProTimeMean[reset], ProTimeStd[reset])
          # EventTime[-1][reset] = Clock + new_times*(new_times>0) #max(0,ProTimeMean[i]+rnorm(1)*ProTimeStd[i])
        }
        
        
      }else{ # stock has arrived
        ComponentType = event-1;
        Inventory[ComponentType] = Inventory[ComponentType]+1
        if (Inventory[ComponentType]>=BaseStockLevel[ComponentType]){
          EventTime[event]=1e5
          if (Clock>WarmUp){ TotalCost=TotalCost+(Clock-OldClock)* sum(OldInventory*HoldingCost) }
        }
      }
    }
    fn = (TotalProfit-TotalCost)/(TotalTime-WarmUp);
    fnSum = fnSum + fn;
  }
  
  fnAvg = fnSum / runlength
  
  return(fnAvg)
}

old_Build_ATO_TestFun = function(seed, numtestseeds=50){
  trainseed = seed*10000 + numtestseeds
  testseeds = seed*10000 + 1:numtestseeds
  
  
  
  TestFun_i = function(xs){
    # the "TruePerf" is when seed=0
    if(xs[9]==0){
      Test_outputs = sapply(testseeds, function(si)RunATO(xs[1:8], seed=si))
      return(mean(Test_outputs))
      
    } else if (xs[9]>0){
      return( RunATO(xs[1:8], seed=trainseed+xs[9]) )
    }
  }
  
  TestFun = function(xs){
    xs = Check_X(xs, 8, T, "ATO TestFun")
    apply(xs, 1, TestFun_i)
  }
  
  return(TestFun)
}

Build_Xie_ATO_Testfun = function(seed=1, numtestseeds=50, runlength=5){
  
  trainseed = seed*10000 + numtestseeds
  testseeds = seed*10000 + 1:numtestseeds
  
  # price, holding cost, avg prod time, std prod time, capacity
  items = c(  1,      2,      .15,   .0225,   20,
              2,      2,      .40,    .06,    20,
              3,      2,      .25,    .0375,  20,
              4,      2,      .15,    .0225,  20,
              5,      2,      .25,    .0375,  20,
              6,      2,      .08,    .012,   20,
              7,      2,      .13,    .0195,  20,
              8,      2,      .40,    .06,    20)
  items = matrix(items, 8, 5, byrow=T)
  
  products = c( 3.6,    1,  0,  0,  1,  0,  1,  1,  0,
                3,      1,  0,  0,  0,  1,  1,  1,  0,
                2.4,    0,  1,  0,  1,  0,  1,  0,  0,
                1.8,    0,  0,  1,  1,  0,  1,  0,  1,
                1.2,    0,  0,  1,  0,  1,  1,  1,  0)
  products = matrix(products, 5, 9, byrow=T)
  
  # Both matrices as defined in problem statement
  
  # Any proposed solution must satisfy bk<=ck
  nItems = nrow(items)
  nProducts=5;
  numberkey=6;
  numbernk=nItems-numberkey;            # Num of products, key items and non key items respectively.
  Tmax=70                                         #Length of simulation
  
  
  
  nGen = 10*Tmax*round(sum(products[,1])) # upper bound on number of generated outputs
  
  Simulate = function(x, seed){
    
    Profit=rep(0, runlength);
    ProfitAfter20= rep(0, runlength);
    set.seed(seed)
    bk = round(x)
    
    # Set the substream to the "seed"
    # ArrivalStream.Substream = seed;
    # ProductionKeyStream.Substream = seed;
    # ProductionNonKeyStream.Substream = seed;
    
    # Generate random data
    # OldStream = RandStream.setGlobalStream(ArrivalStream); # Temporarily store old stream
    
    # Arrival=zeros(nProducts,nGen,runlength);
    Arrival = array(0, dim=c(nProducts, nGen, runlength))
    #Generate time of next order arrival per product
    for (k in 1:nProducts) Arrival[k,,]= -log(runif(nGen*runlength))*(1/products[k,1]) # exprnd(1/products[k,1],1,nGen,runlength);
    
    
    
    # Generate production times
    # RandStream.setGlobalStream(ProductionKeyStream);
    ProdTimeKey= matrix(rnorm(nGen*runlength), nGen, runlength) #normrnd(0,1,nGen,runlength);
    
    # Generate Uniforms used to generate orders
    # RandStream.setGlobalStream(ProductionNonKeyStream);
    ProdTimeNonKey= matrix(rnorm(nGen*runlength), nGen, runlength) #normrnd(0,1,nGen,runlength);
    
    # Restore old random number stream
    # RandStream.setGlobalStream(OldStream);
    
    for (k in 1:runlength){
      # Initialize this replication
      Inventory = bk;   # Tracks available inventory for each item
      Orders    = Arrival[,1,k] # Next order arrival times
      
      
      # ItemAvailability contains the time at which a replenishment order for a
      # given item(row) will be ready. Each order replenishes one unit of
      # one item.
      
      itemAvail = lapply(1:nItems, function(a)0)
      A = rep(2, nProducts) # indeces ticking through each row of Arrival (times)
      prevMinTime=0      #Holds the time of previous order fulfilled (to calculate holding cost)
      p=1 # Index for Key item production times
      q=1 # Index for non-Key production times
      
      # Main simulation:
      #**** Loop through orders, as they happened (smallest time first) and identify ****
      #**** whether key and non-key items are available, calculate profit, etc.      ****
      
      # While there are orders to be satisfied.
      while(min(Orders)<=Tmax){
        # find next order to satisfy and when it happens
        minTime = min(Orders)
        minProd = which.min(Orders)
        
        # generate time of next order and iterate the Arrival time 'A' index
        Orders[minProd] = Orders[minProd] + Arrival[minProd, A[minProd], k]
        A[minProd] = A[minProd] + 1
        
        # Add inventory that have become available upto minTime
        Inventory = Inventory + sapply(itemAvail, function(ir)sum(0<ir&ir<=minTime))
        itemAvail = lapply(itemAvail, function(ir)ir[ir>minTime])
        
        # if(all key items available) make product, etc
        keyavail = all(products[minProd, 1+1:numberkey]<=Inventory[1:numberkey])
        
        
        
        ##################################################################################################
        ##################################################################################################
        ##################################################################################################
        
        
        if(keyavail){ # all key items available
          
          Profit[k] = Profit[k] + sum(products[minProd,1+1:numberkey]*items[1:numberkey,1])
          
          if(minTime>=20)ProfitAfter20[k] = ProfitAfter20[k] + sum(products[minProd,1+1:numberkey]*items[1:numberkey,1]);
          
          for (r in 1:numberkey){
            # Decrease inventory and place replenishment orders for the amount of key items used
            num = products[minProd,r+1]
            if (num!=0){
              Inventory[r]=Inventory[r]-num;
              for (g in 1:num){
                tt = length(itemAvail[[r]])
                itemAvail[[r]][tt+1] = max(minTime,itemAvail[[r]][tt]) +(items[r,3] + items[r,4]*ProdTimeKey[p,k])
                p=p+1
              }
            }
          }
          
          # For each non-key item available, use it, decrease inventory, increase profit and place replenishment order.
          for (j in 1:numbernk){
            r   = numberkey +j
            num = products[minProd,r+1]
            
            if(num<=Inventory[r] && num!=0){
              
              Profit[k]=Profit[k]+items[r,1]*num;
              if(minTime>=20) ProfitAfter20[k]=ProfitAfter20[k]+items[r,1]*num
              
              Inventory[r] = Inventory[r]-num
              for (g in 1:num){
                tt = length(itemAvail[[r]])
                itemAvail[[r]][tt+1]=max(minTime,itemAvail[[r]][tt])+(items[j+numberkey,4]*ProdTimeNonKey[q,k]+items[j+numberkey,3])
                q=q+1;
              }
            }
          }
          
          # itemAvailability updated and resized!
        }
        
        
        ##################################################################################################
        ##################################################################################################
        ##################################################################################################
        
        
        
        Profit[k]=Profit[k]-sum(Inventory*items[,2])*(minTime-prevMinTime);
        if(minTime>=20) ProfitAfter20[k] = ProfitAfter20[k]-sum(Inventory*items[,2])*(minTime-prevMinTime);
        
        prevMinTime=minTime;
      }
    }
    
    dailyProfitAfter20 = ProfitAfter20/(Tmax-20);
    fn    = mean(dailyProfitAfter20);
    FnVar = var(dailyProfitAfter20)/runlength;
    # print(fn)
    
    fn
  }
  
  TestFun_i = function(xs){
    # the "TruePerf" is when seed=0
    if(any(xs[1:8]<0) | any(20<xs[1:8]))stop("Xie ATO x is out of bounds") 
    
    if(xs[9]==0){
      Test_outputs = sapply(testseeds, function(si)Simulate(x=xs[1:8], seed=si))
      return(mean(Test_outputs))
      
    } else if (xs[9]>0){
      return( Simulate(xs[1:8], seed=trainseed+xs[9]) )
    }
  }
  
  TestFun = function(xs){
    xs = Check_X(xs, 8, T, "Xie ATO TestFun")
    apply(xs, 1, TestFun_i)
  }
  
  return(TestFun)
  
}

Build_Xie_ATO_cpp_Testfun_old = function(baseseed=1, numtestseeds=200, runlength=5){
  Rcpp::sourceCpp('TestFuns/ATO.cpp')
  items = c(  1,      2,      .15,   .0225,   20,
              2,      2,      .40,    .06,    20,
              3,      2,      .25,    .0375,  20,
              4,      2,      .15,    .0225,  20,
              5,      2,      .25,    .0375,  20,
              6,      2,      .08,    .012,   20,
              7,      2,      .13,    .0195,  20,
              8,      2,      .40,    .06,    20)
  items = matrix(items, 8, 5, byrow=T)
  
  products = c( 3.6,    1,  0,  0,  1,  0,  1,  1,  0,
                3,      1,  0,  0,  0,  1,  1,  1,  0,
                2.4,    0,  1,  0,  1,  0,  1,  0,  0,
                1.8,    0,  0,  1,  1,  0,  1,  0,  1,
                1.2,    0,  0,  1,  0,  1,  1,  1,  0)
  products  = matrix(products, 5, 9, byrow=T)
  Tmax      = 70
  nProducts = 5
  nGen      = 10*Tmax*round(sum(products[,1]))
  
  Make_Stream = function(s, runlength, baseseed){
    cat("calling makestream\n")
    set.seed(1000*baseseed+s)
    Arrival = array(0, dim=c(runlength, nProducts, nGen))
    for (k in 1:nProducts) Arrival[,k,]= -log(runif(nGen*runlength))*(1/products[k,1]) 
    
    ProdTimeKey= matrix(rnorm(nGen*runlength), runlength, nGen)
    ProdTimeNonKey= matrix(rnorm(nGen*runlength), runlength, nGen)
    
    list(Arrival=Arrival, ProdTimeKey=ProdTimeKey, ProdTimeNonKey=ProdTimeNonKey, runlength=runlength)
  }
  
  RV = lapply(1:10001, function(a)NULL)
  RV[[1]] = Make_Stream(0, numtestseeds*runlength, baseseed)
  
  simcalls = 0
  
  TestFun_i = function(xs){
    
    ss = xs[9]+1
    # Generate new RNG streams
    # if(length(RV)< ss){
    #   RV[[ss]] <<- Make_Stream(ss, runlength, baseseed)
    # }else 
    if(is.null(RV[[ss]])){
      cat(ss, is.null(RV[[ss]]), "\n")
      RV[[ss]] <<- Make_Stream(ss, runlength, baseseed)
      # RV <<- RV
    }
    
    out = sapply(1:RV[[ss]]$runlength, function(k){
      simcalls <<- simcalls+1
      Simulate_ATO(
        xs[1:8], 
        RV[[ss]]$Arrival[k,,], 
        RV[[ss]]$ProdTimeKey[k,], 
        RV[[ss]]$ProdTimeNonKey[k,], 
        items, 
        products)})
    return(mean(out))
  }
  
  TestFun = function(xs){
    xs = Check_X(xs, 8, T, "ATO cpp")
    out = 1:nrow(xs)
    for(i in 1:nrow(xs))out[i] = TestFun_i(xs[i,])
    # cat("\n",length(RV))
    out
  }
  
  Get_RV = function() RV
  
  Get_simcalls = function()simcalls
  
  c(TestFun, Get_RV, Get_simcalls)
}

Build_Xie_ATO_cpp_Testfun = function(baseseed=1, numtestseeds=500, runlength=5){
  Rcpp::sourceCpp('TestFuns/ATO.cpp')
  items = c(  1,      2,      .15,   .0225,   20,
              2,      2,      .40,    .06,    20,
              3,      2,      .25,    .0375,  20,
              4,      2,      .15,    .0225,  20,
              5,      2,      .25,    .0375,  20,
              6,      2,      .08,    .012,   20,
              7,      2,      .13,    .0195,  20,
              8,      2,      .40,    .06,    20)
  items = matrix(items, 8, 5, byrow=T)
  
  products = c( 3.6,    1,  0,  0,  1,  0,  1,  1,  0,
                3,      1,  0,  0,  0,  1,  1,  1,  0,
                2.4,    0,  1,  0,  1,  0,  1,  0,  0,
                1.8,    0,  0,  1,  1,  0,  1,  0,  1,
                1.2,    0,  0,  1,  0,  1,  1,  1,  0)
  products  = matrix(products, 5, 9, byrow=T)
  Tmax      = 70
  nProducts = 5
  nGenA     = 5*Tmax*round(max(products[,1]))
  nGen      = 5*Tmax*round(sum(products[,1]))
  
  Make_Stream = function(s, runlength, baseseed){
    # cat("calling makestream\n")
    set.seed(1000*baseseed+s)
    Arrival = lapply(1:runlength, function(i){
      t(sapply(1:nProducts, function(k) -log(runif(nGenA))*(1/products[k,1])))
    })
    ProdTimeKey= lapply(1:runlength, function(i)rnorm(nGen))
    ProdTimeNonKey= lapply(1:runlength, function(i)rnorm(nGen))
    
    list(Arrival=Arrival, ProdTimeKey=ProdTimeKey, ProdTimeNonKey=ProdTimeNonKey, runlength=runlength)
  }
  
  # RV = lapply(1:10001, function(a)NULL)
  TestStreams = Make_Stream(0, numtestseeds*runlength, baseseed)
  
  simcalls = 0
  
  TestFun_i = function(xs){
    
    if(xs[9]==0){
      RV = TestStreams
    }else{
      # Generate new RNG streams
      RV = Make_Stream(xs[9], runlength, baseseed)
    }
    
    out = sapply(1:RV$runlength, function(k){
      simcalls <<- simcalls+1
      Simulate_ATO(
        xs[1:8], 
        RV$Arrival[[k]], 
        RV$ProdTimeKey[[k]], 
        RV$ProdTimeNonKey[[k]], 
        items, 
        products)})
    return(mean(out))
  }
  
  TestFun = function(xs){
    cat("\nRunning ATO sim...")
    if (is.null(dim(xs)))xs=matrix(xs,1)
    xs = Check_X(xs, 8, T, "ATO cpp")
    out = 1:nrow(xs)
    for(i in 1:nrow(xs))out[i] = TestFun_i(xs[i,])
    cat("done, ")
    out
  }
  
  Get_RV = function() TestStreams
  
  Get_simcalls = function()simcalls
  
  attr(TestFun, 'ran') = rbind(rep(0,8), rep(20,8))
  attr(TestFun, 'name') = "ATO"
  
  
  c(TestFun, Get_RV, Get_simcalls)
  # TestFun
}

Build_Matlab_ATO = function(BOseed, numtestseeds=100, runlength1=5, HN1=F){
  
  # stop("gotta do some sanity checking to see this works as expected!")
  BaseSeed = 10000*BOseed
  
  Matlab_ATO = function(x, runlength, path_to_ATO_batch="/home/maths/phrnaj/ATO", HN=HN1){
    
    folder = getwd()
    if(is.null(dim(x))&length(x)==9) x = matrix(x,1)
    
    out_fname = tempfile(fileext = ".txt")
    mat_fname = tempfile(fileext = ".m")
    
    row_str = sapply(1:nrow(x), function(a)paste(x[a,],collapse=","))
    mat_str = paste(row_str, collapse="; ")
    
    runl_str = paste(runlength, collapse=",")
    
    if(HN){
      code = paste( "cd('",path_to_ATO_batch,"'); \np=ATOHongNelson_batch([", mat_str, "],[",runl_str,"]); \nsave('",out_fname,"','p','-ascii')", sep="")
    }else{
      code = paste( "cd('",path_to_ATO_batch,"'); \np=ATO_batch([", mat_str, "],[",runl_str,"]); \nsave('",out_fname,"','p','-ascii')", sep="")
    }
    write(code, mat_fname)
    
    exec1 = "matlab  -nojvm -nodesktop -nosplash -nodisplay -r  \"try, run('"
    exec2 = "'); catch err, disp(err.message); exit(1); end; exit(0);\""
    full_exec = paste(exec1, mat_fname, exec2, sep="")
    
    cat("calling MATLAB...")
    system(full_exec, wait=T, ignore.stdout=T)
    cat("done, ")
    as.numeric(readLines(out_fname))
  }
  
  TestFun = function(x){
    x = Check_X(x,8,T,"matlab ATO")
    
    runl = (1 + (x[,9]==0)*(numtestseeds-1)) * runlength1
    
    x[,9] = x[,9] + BaseSeed
    
    output = Matlab_ATO(x, runl)
    
    return(output)
    
  }
  
  return(TestFun)
  
}

Build_ATO_TestFun = function(seed, numtestseeds=50, runlength=5){
  
  trainseed = seed*10000 + numtestseeds
  testseeds = seed*10000 + 1:numtestseeds
  
  
  # Args:
  #   BaseStockLEvel: target stock quantities for each product
  #   runlength: integer, repetitions of the simulation
  #   seed: integer, to force deterministic simulation
  #
  # Returns:
  #   fnAvg: Average profit over runlength reps of simulation
  
  NumComponentType=8;               #% number of component types
  
  ProTimeMean = c(0.15, 0.4, 0.25, 0.15, 0.25, 0.08, 0.13, 0.4)
  ProTimeStd = 0.15*ProTimeMean
  Profit = 1:8
  HoldingCost = matrix(2, 1, NumComponentType)
  # % holding cost of each component in
  # % inventory
  # 
  # % Parameters of Customers
  
  ArrivalRate = 12             #% assume Possion arrival
  iArrivalRate = 1/ArrivalRate
  # %NumCustomerType=5;             % number of customer types
  # %CustomerProb=[0.3,0.25,0.2,0.15,0.1];   % probability of each customer
  KeyComponent=c(1,0,0,1,0,1,0,0,
                 1,0,0,0,1,1,0,0,
                 0,1,0,1,0,1,0,0,
                 0,0,1,1,0,1,0,0,
                 0,0,1,0,1,1,0,0)
  KeyComponent = matrix(KeyComponent, nrow=5, byrow=TRUE)
  
  NonKeyComponent=c(0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,1,
                    0,0,0,0,0,0,1,0)
  
  NonKeyComponent = matrix(NonKeyComponent, nrow=5, byrow=TRUE)
  
  WarmUp = 20
  TotalTime = 70
  
  Simulate = function(BaseStockLevel=rep(1,8), seed=1){
    # BaseStockLevel = floor(BaseStockLevel)
    #% Simulation
    fnSum = 0.0
    set.seed(seed)
    Customer_series = sample(size=TotalTime*ArrivalRate*runlength*4, 1:5, prob=c(0.3,0.25,0.2,0.15,0.1), replace=T)
    t = 1
    
    for (k in 1:runlength){
      
      # runif_preset = runif()
      
      EventTime = matrix(1e5, 1, 1+NumComponentType)
      EventTime[1] = -log(runif(1))*iArrivalRate
      TotalProfit = 0
      TotalCost = 0
      Inventory = ceiling(BaseStockLevel)
      
      Clock = 0
      
      while (Clock<TotalTime){
        
        OldInventory = Inventory
        OldClock = Clock
        Clock = min(EventTime)
        event = which.min(EventTime)
        if (event==1){ # a customer has come!
          # Reset clock
          EventTime[1]=Clock-log(runif(1))*iArrivalRate
          
          # Get new customer
          CustomerType = Customer_series[t]
          t=t+1
          
          KeyOrder = KeyComponent[CustomerType,]
          NonKeyOrder = NonKeyComponent[CustomerType,]
          
          Sell = all(Inventory>=KeyOrder)
          NonKeyOrder = Inventory * (Inventory<=NonKeyOrder) + NonKeyOrder*(Inventory>NonKeyOrder)
          
          if (Sell){
            Inventory = Inventory - KeyOrder - NonKeyOrder
            if (Clock>WarmUp){ TotalProfit = TotalProfit + sum(Profit*(KeyOrder + NonKeyOrder)) }
            
            for (i in 1:NumComponentType){
              if (Inventory[i]<BaseStockLevel[i]  && EventTime[i+1]>1e4){
                EventTime[i+1] = Clock + max(0,ProTimeMean[i]+rnorm(1)*ProTimeStd[i])
              }
            }
            # reset = ((Inventory<BaseStockLevel)  & (EventTime[-1]>1e4))
            # new_times = rnorm(sum(reset), ProTimeMean[reset], ProTimeStd[reset])
            # EventTime[-1][reset] = Clock + new_times*(new_times>0) #max(0,ProTimeMean[i]+rnorm(1)*ProTimeStd[i])
          }
          
          
        }else{ # stock has arrived
          ComponentType = event-1;
          Inventory[ComponentType] = Inventory[ComponentType]+1
          if (Inventory[ComponentType]>=BaseStockLevel[ComponentType]){
            EventTime[event]=1e5
            if (Clock>WarmUp){ TotalCost=TotalCost+(Clock-OldClock)* sum(OldInventory*HoldingCost) }
          }
        }
      }
      fn = (TotalProfit-TotalCost)/(TotalTime-WarmUp);
      fnSum = fnSum + fn;
    }
    
    fnAvg = fnSum / runlength
    
    return(fnAvg)
  }
  
  TestFun_i = function(xs){
    # the "TruePerf" is when seed=0
    if(any(xs[1:8]<0) | any(20<xs[1:8]))stop("ATO x is out of bounds") 
    
    if(xs[9]==0){
      Test_outputs = sapply(testseeds, function(si)Simulate(xs[1:8], seed=si))
      return(mean(Test_outputs))
      
    } else if (xs[9]>0){
      return( Simulate(xs[1:8], seed=trainseed+xs[9]) )
    }
  }
  
  TestFun = function(xs){
    xs = Check_X(xs, 8, T, "ATO TestFun")
    apply(xs, 1, TestFun_i)
  }
  
  # return(TestFun)
  return(TestFun)
  
}

Build_GP_TestFun = function(seed, dims, TPars){
  # Args:
  #   seed: the fixed random number generator seed
  #   dims: dimension of decision variable
  #   TPars: 5 hyperparams, dims*L_theta, S_theta, dims*L_eps, S_eps, C_eps
  #
  # Returns:
  #   TestFun: a function that takes continuous input and a seed, returns scalar
  
  ThetaFuns    = GenFun(LX=rep(TPars[1],dims), SS2=TPars[2], SC2=0, seed=1000*seed)
  Theta        = ThetaFuns$fun
  EPS          = lapply(1:500, function(i)GenFun(LX = rep(TPars[3], dims), SS2 = TPars[4], SC2=TPars[5], seed = seed*1000+i)$fun)
  
  DD = 1:dims
  SS = dims+1
  
  TestFun_i = function(xsi) Theta(xsi[DD]) + if(xsi[SS]>0&xsi[SS]<=length(EPS)) EPS[[xsi[SS]]](xsi[DD]) else 0
  
  TestFun  = function(xs){
    xs = Check_X(xs, dims, T, "GP TestFun")
    apply(xs, 1, TestFun_i)
  }
  
  return(list(fun=TestFun, grad=ThetaFuns$dfun))
}

########################################################################################################
# Ambulance

Build_Ambulance_Testfun = function(baseseed=1, numtestseeds=10000, runlength=5, NumCalls = 30){
  Rcpp::sourceCpp('TestFuns/Ambulances.cpp')
  cat(" Ambulances.cpp compilation complete. ")
  
  simcalls = 0 
  Make_Stream = function(s, runlength, baseseed){
    # cat("calling makestream\n")
    set.seed(10000*baseseed+s)
    
    # NumCalls = 
    
    CallTimes = lapply(1:runlength, function(r)cumsum( rexp(NumCalls, rate = 1/60) ))
    
    CallLocs  = lapply(1:runlength, function(r){
      CallLocs_r = matrix(0, 0, 2)
      while (nrow(CallLocs_r)<NumCalls){
        Calls = round(2*( NumCalls - nrow(CallLocs_r) ))
        u = matrix(runif(3*Calls), Calls, 3)
        Acc = 1.6*u[,3] <= 1.6 - abs(u[,1]-0.8) - abs(u[,2]-0.8)
        CallLocs_r = rbind(CallLocs_r, u[Acc, 1:2])
      }
      CallLocs_r = CallLocs_r[1:NumCalls,]
    })
    
    Stimes    = lapply(1:runlength, function(r)rgamma(NumCalls, shape=9, scale = 1/12))
    
    list(CallTimes=CallTimes, CallLocs=CallLocs, Stimes=Stimes, runlength=runlength)
    
  }
  
  TestStreams = Make_Stream(0, numtestseeds*runlength, baseseed)
  
  TestFun_i = function(xs){
    
    if(xs[7]==0){
      RV = TestStreams
    }else{
      # Generate new RNG streams
      RV = Make_Stream(xs[7], runlength, baseseed)
    }
    
    xs = matrix(xs[1:6], 3,2)*0.05
    
    simcalls <<- simcalls + RV$runlength
    
    out = sapply(1:RV$runlength, function(k){
      Ambulances_Square(xs, RV$CallTimes[[k]], RV$CallLocs[[k]], RV$Stimes[[k]])
    })
    return(-mean(out))
  }
  
  TestFun = function(xs){
    cat("\nRunning Ambulance sim...")
    if (is.null(dim(xs)))xs = matrix(xs, 1)
    xs = Check_X(xs, 6, T, "Ambulance cpp")
    out = 1:nrow(xs)
    for(i in 1:nrow(xs))out[i] = TestFun_i(xs[i,])
    cat("done, ")
    out
  }
  
  Get_RV = function() TestStreams
  
  Get_simcalls = function()simcalls
  
  attr(TestFun, 'ran') = rbind(rep(0,6), rep(20,6))
  attr(TestFun, 'name') = "Ambulances"
  
  c(TestFun, Get_RV, Get_simcalls)
  # TestFun
}

Build_DailyAmbulance_Testfun = function(baseseed=1, numtestseeds=10000, runlength=5, Simtime=1800, reduce_out=mean){
  Rcpp::sourceCpp('TestFuns/Ambulances.cpp')
  cat(" Ambulances.cpp compilation complete. ")
  
  simcalls = 0 
  Make_Stream = function(s, runlength, baseseed){
    # cat("calling makestream\n")
    set.seed(10000*baseseed+s)
    
    
    NumCalls = 2 * 30 * Simtime/1800
    
    CallTimes = lapply(1:runlength, function(r)cumsum( rexp(NumCalls, rate = 1/60) ))
    
    CallLocs  = lapply(1:runlength, function(r){
      CallLocs_r = matrix(0, 0, 2)
      while (nrow(CallLocs_r)<NumCalls){
        Calls = round(2*( NumCalls - nrow(CallLocs_r) ))
        u = matrix(runif(3*Calls), Calls, 3)
        Acc = 1.6*u[,3] <= 1.6 - abs(u[,1]-0.8) - abs(u[,2]-0.8)
        CallLocs_r = rbind(CallLocs_r, u[Acc, 1:2])
      }
      CallLocs_r = CallLocs_r[1:NumCalls,]
    })
    
    Stimes    = lapply(1:runlength, function(r)rgamma(NumCalls, shape=9, scale = 1/12))
    
    list(CallTimes=CallTimes, CallLocs=CallLocs, Stimes=Stimes, runlength=runlength)
    
  }
  
  TestStreams = Make_Stream(0, numtestseeds*runlength, baseseed)
  
  TestFun_i = function(xs){
    
    if(xs[7]==0){
      RV = TestStreams
    }else{
      # Generate new RNG streams
      RV = Make_Stream(xs[7], runlength, baseseed)
    }
    
    xs = matrix(xs[1:6], 3,2)*0.05
    
    simcalls <<- simcalls + RV$runlength
    
    out = sapply(1:RV$runlength, function(k){
      O = Ambulances_Square_timed(xs, RV$CallTimes[[k]], RV$CallLocs[[k]], RV$Stimes[[k]], Simtime)
      sO = sort(O[O>0])
      
      quantile = sO[round(0.75*length(sO))]
      return(quantile)
    })
    return(-mean(out))
  }
  
  TestFun = function(xs){
    cat("\nRunning Ambulance sim...")
    if (is.null(dim(xs)))xs = matrix(xs, 1)
    xs = Check_X(xs, 6, T, "Ambulance cpp")
    out = 1:nrow(xs)
    for(i in 1:nrow(xs))out[i] = TestFun_i(xs[i,])
    cat("done, ")
    out
  }
  
  Get_RV = function() TestStreams
  
  Get_simcalls = function()simcalls
  
  attr(TestFun, 'ran') = rbind(rep(0,6), rep(20,6))
  attr(TestFun, 'name') = "Ambulances"
  
  c(TestFun, Get_RV, Get_simcalls)
  # TestFun
}

########################################################################################################
# GP

GenFun = function(seed=NULL, LX=20, SS2=25, SC2=100, Nx=NULL){
  # Create a single continuous GP test function  
  # Testfun:[0,100]^dims * N^+ -> R where the seed=0
  # gives the ground truth.
  # 
  # Args
  #   seed: fixed seed
  #   LX: SE kernel length scales, this determines dimensions
  #   SS2: the signal variance of the function
  #   SC2: the variance of global seed offsets
  # 
  # Returns
  #   result: a continouos callable function: R^d -> R
  
  if(is.null(Nx)){Nx = prod(200/LX); Nx = min(Nx, 1000)}
  
  dims = length(LX)
  ilX  = -0.5/LX^2
  SEKernel = function(x1,x2){
    x1 = matrix(x1,ncol=dims)
    x2 = matrix(x2,ncol=dims)
    exp(Reduce(f = "+",lapply(1:dims,function(d)outer(x1[,d],x2[,d],"-")^2*ilX[d]  )))
  }
  
  set.seed(seed)
  
  ran    = matrix(c(0, 100), 2, dims)
  XX0    = LHCran(Nx, ran)
  TrueK  = SEKernel(XX0,XX0) #+ 0.0001*diag(nrow(XX0))
  TrueiK = rcppeigen_invert_matrix(TrueK)
  
  SC0  = rnorm(1,0,sqrt(SC2))
  y    = mvrnorm(1,rep(0,nrow(XX0)), TrueK)
  iKy  = (TrueiK%*%y)[,1]
  SS1  = sqrt(SS2)
  
  result=function(X){
    
    if (is.null(dim(X)) & length(X)==dims) X=matrix(X,ncol=dims, nrow=1)
    
    if (length(X)%%dims!=0)stop("wrong shape input to testfun")
    
    if (ncol(X)!=dims)stop("wrong shape input to GP testfun")
    
    Ks = SEKernel(X,XX0)
    SC0 + SS1*as.numeric(Ks%*%iKy)
  }
  
  dresult = function(X){
    X  = Check_X(X, dims, F, "dtestfun")
    Ks = SEKernel(X,XX0)
    2*sapply(1:dims, function(d)(outer(X[,d], XX0[,d], "-")*ilX[d]*Ks )%*%iKy)*SS1
  }
  
  attr(result, 'ran') = ran
  
  return(list(fun=result, dfun=dresult))
}

Build_Discrete_GP_CompSph = function(seed, dims, TPars){
  # Args:
  #   seed: the fixed random number generator seed
  #   dims: dimension of decision variable
  #   TPars: 4 hyperparams, dims*L_theta, S_theta, S_eps, C_eps
  #
  # Returns:
  #   TestFun: a function that takes continuous input and a seed, returns scalar
  
  ThetaFuns    = GenFun(LX=rep(TPars[1],dims), SS2=TPars[2], SC2=0, seed=1000*seed)
  Theta        = ThetaFuns$fun
  offsets      = c(0, rnorm(1000)) * sqrt(TPars[length(TPars)])
  noises       = c(0, rep(1, 1000)) * sqrt(TPars[length(TPars)-1])
  
  DD = 1:dims
  SS = dims+1
  
  TestFun_i = function(xsi) {
    s = xsi[length(xsi)]
    s = (s-1)%%1000 + 1
    xsi[length(xsi)] = s
    Theta(xsi[DD]) + rnorm(1, offsets[xsi[SS]+1], noises[xsi[SS]+1])
  }
  
  TestFun  = function(xs){
    xs = Check_X(xs, dims, T, "GP Disc TestFun")
    apply(xs, 1, TestFun_i)
  }
  
  attr(TestFun, 'ran') = rbind(rep(0, dims), rep(100, dims))
  
  return(list(fun=TestFun, grad=ThetaFuns$dfun))
}

Build_Discrete_GP_CompSph2 = function(seed, TPars){
  # Args:
  #   seed: the fixed random number generator seed
  #   TPars: 4 hyperparams, L_theta, S_theta, S_eps, C_eps
  #
  # Returns:
  #   TestFun: a function that takes continuous input and a seed, returns scalar
  
  LX = TPars[1:(length(TPars)-3)]
  dims = length(LX)
  SS2 = TPars[length(TPars)-2]
  
  ThetaFuns    = GenFun(LX=LX, SS2=SS2, SC2=0, seed=1000*seed)
  Theta        = ThetaFuns$fun
  noises       = c(0, rep(1, 10000)) * sqrt(TPars[length(TPars)-1])
  offsets      = c(0, rnorm(10000)) * sqrt(TPars[length(TPars)])
  
  DD = 1:dims
  SS = dims+1
  
  TestFun_i = function(xsi) {
    s = xsi[length(xsi)]
    Theta(xsi[DD]) + rnorm(1, offsets[s+1], noises[s+1])
  }
  
  TestFun  = function(xs){
    xs = Check_X(xs, dims, T, "GP Disc TestFun")
    apply(xs, 1, TestFun_i)
  }
  
  return(list(fun=TestFun, grad=ThetaFuns$dfun))
}

Build_Wiggle = function(seed=1, dims=1, LX=10, SS2=100^2, WX=5, WW2=30^2, rho=0.5){
  
  theta = GenFun(seed=seed*1000,
                 LX=rep(LX,dims),
                 SS2=SS2, 
                 SC2=0
                 )$fun
  # browser()
  wiggles = lapply(1:1000, function(i)GenFun(seed = i+seed*1000,
                                             LX=rep(WX, dims),
                                             SS2=(WW2*rho),
                                             SC2=0
                                             )$fun)
  
  # browser()
  noiseSD = sqrt( WW2*(1-rho) )
  
  TestFun_i = function(xs){
    se = xs[dims+1]
    xx = xs[1:dims]
    eps=0
    
    if(se>0){ se=(se-1)%%1000 +1;  eps = wiggles[[se]](xx) + rnorm(1, 0, noiseSD)}
    
    return(theta(xx)+eps)
  }
  
  TestFun = function(XS){
    XS = Check_X(XS, dims, T, "GPwiggle")
    return( apply(XS, 1, TestFun_i) )
  }
  
  attr(TestFun, 'ran') = rbind(rep(0, dims), rep(100, dims))
  
  return(TestFun)
}
