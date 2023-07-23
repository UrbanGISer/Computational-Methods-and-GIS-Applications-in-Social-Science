# CMGIS_Rtool: Computational Methods and GIS Applications in Social Sciences,Third Edition , CRC Press
# coding contributor: Lingbo Liu
# Center for Geographic Analysis, Harvard University
# School of Urban Design, Wuhan University 

tool_exec<- function(in_params, out_params){
  
  #####################################################################################################  
  ### Check/Load Required Packages  
  #####################################################################################################   
  arc.progress_label("Loading packages...")
  arc.progress_pos(0)

  if(!requireNamespace("quadprog", quietly = TRUE))
    install.packages("quadprog", quiet = TRUE) 
  if(!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", quiet = TRUE) 

  require(quadprog)
  require(dplyr)
  
  arc.progress_label("Packages Load Successfully")
  arc.progress_pos(10)
  
  ##################################################################################################### 
  ### Define input/output parameters
  #####################################################################################################
  # OD data
  input_OD_data <- in_params[[1]]
  OID<-in_params[[2]]
  DID<-in_params[[3]]
  Distance<-in_params[[4]]
  # Demand data
  input_demand_data<- in_params[[5]]
  DemandID<-in_params[[6]]
  Demandfield<-in_params[[7]]
  
  # Hospital data
  input_fixed_hosp<- in_params[[8]]
  FixCapcacity<-in_params[[9]]   # field from input_fixed_hosp
  NewCapacity <-in_params[[10]]  # numeric
  
  # Accessibility Mode
  AccModel<-in_params[[11]]      # 2SFCA Exponential Gravity
  # Accessibility Mode Parameter
  Threshold<-in_params[[12]]     # 2SFCA
  Gravitybeta<-in_params[[13]]   # Gravity
  Expon<-in_params[[14]]         # Exponential
  # Output table  
  Out_supply_data <- out_params[[1]]
  Out_Acc_data  <- out_params[[2]]
  
  ##################################################################################################### 
  ### Load Data and Create Dataframe R Object 
  #####################################################################################################
  arc.progress_label("Loading data...")
  arc.progress_pos(30)

  # Read Data
  ODdata <- arc.open(input_OD_data)
  demanddata <- arc.open(input_demand_data)
  hospdata <- arc.open(input_fixed_hosp)

  # Extract List
  OD_list<-c(OID,DID,Distance)
  Demand_list<-c(DemandID,Demandfield)

  # update data
  ODdata <- arc.select(ODdata,fields = OD_list)
  demanddata <-arc.select(demanddata,fields = Demand_list)
  names(ODdata) <-c("OID","DID","Dist")
  names(demanddata) <-c("OID","Popu")
  hospdata <- arc.open( input_fixed_hosp)
  hosp <- arc.select(hospdata)
  fixhosp<-subset(hosp,select=FixCapcacity)
  names(fixhosp)<-c('Supply')
  # Get important parameters, row and column
  
  demandpt <-nrow(demanddata )
  supplypt  <- nrow(ODdata)/demandpt
  demand_popu <- sum(demanddata$Popu)
  fixH<-nrow(fixhosp)
  fixcapacity<-sum(fixhosp$Supply)
  Totalcapacity<-fixcapacity+NewCapacity
  newH<-supplypt-fixH
  ave_accessibility <- Totalcapacity / demand_popu 
  
  print(paste(demandpt,'demand locations with total population of ',demand_popu))
  print(paste(supplypt,'facilites with total capacity of',round(Totalcapacity,3)))
  print(paste( fixH,'fixed facilities with total capacity of',round(fixcapacity,3)))
  print(paste( newH, 'new facilities  with total capacity of',NewCapacity)) 
  print(paste( 'Average Accessibility Score is',ave_accessibility))  
  
  ##################################################################################################### 
  ### Prepare Matrix data for subsequent analysis
  #####################################################################################################
  arc.progress_label("Building Matrix...")
  arc.progress_pos(50) 
  
  
  ODpotent<-merge(ODdata,demanddata,by='OID')
  
  # Calculate distance decay effect f(dij) in 3 modes
  if (AccModel== '2SFCA') {
    ODpotent$fdij<-ifelse(ODpotent$Dist<= Threshold, 1, 0)
  } else if (AccModel== 'Gravity') {
    ODpotent$fdij<-ODpotent$Dist^(-1*Gravitybeta)
  } else {
    ODpotent$fdij<-exp(-1*ODpotent$Dist*Expon)
  }
  
  # Calculate distance-decay-effect weighted population D*(f(dij))
  ODpotent$Dfdkj<-ODpotent$Popu*ODpotent$fdij
  
  # Calculate total distance-weighted demand population of each facility
  Sum_Dfdki<-ODpotent %>% group_by(DID)%>%summarise(Sum_Dfdki=sum(Dfdkj))
  
  # Calculate Fij matrix 
  ODpotent<-merge(ODpotent,Sum_Dfdki,by='DID')
  ODpotent$Fij<-ODpotent$fdij/ODpotent$Sum_Dfdki
  
  #----------------------------------------------------------#
  # Prepare Matrix for Qp package
  ODpotent<-ODpotent[order(ODpotent[,'OID'],ODpotent[,'DID']),]
  vec_Fij<-ODpotent[,"Fij"]
  vec_D<-demanddata[,"Popu"]
  
  # Prepare F D for H=F'DF
  Fij<-t(matrix(vec_Fij,nrow=supplypt,ncol=demandpt))
  D<-matrix(diag(vec_D),ncol=demandpt) 
  
  #------Define ones and zeros function in matlab-----------#
  newone<-function(m,n) {
    if (m==1) return (matrix(rep(1, n),ncol = n))
    else if (n==1) return (matrix(rep(1, m),ncol = 1))
  }
  
  newzero<-function(m,n) {
    if (m==1) return (matrix(rep(0, n),ncol = n))
    else if (n==1) return (matrix(rep(0, m),ncol = 1))
  }
  
  # Average Accessibility
  A <-newone(demandpt,1) %*% ave_accessibility
  
  # H=F'DF  f=F'DA
  Dmat = t(Fij)%*% D %*% Fij
  dvec = t(Fij) %*% D %*% A
  
  
  # Ex=d ,Cx<=b, here, B=E, x=S, d=supply_popu
  B <- newone(1,supplypt)
  AA <-diag(c(newone(supplypt,1)))
  lb <- newzero(supplypt,1)
  
  # build  constraints
  Amat <- t(rbind(B,AA))
  bvec <- c(Totalcapacity,lb)
  


  ##################################################################################################### 
  ### Run Rglpk package
  #####################################################################################################
  arc.progress_label("Apply Rglpk...")
  arc.progress_pos(60)
  
  # use all hospitals 
  qp <- solve.QP(Dmat, dvec, Amat, bvec=bvec, meq = 1)
  
  # use fixed hospitals 
  AAfix <-AA[1:fixH,]
  AAsolve <- matrix(c(rep (0 , fixH ) , rep (1 , newH ) ),nrow = 1)
  bvec1  <- c(fixhosp$Supply,NewCapacity)
  Amat1  <- t(rbind(AAfix,AAsolve))
  qp1 <- solve.QP(Dmat, dvec, Amat1, bvec=bvec1, meq = 1)
  
  # Gather data
  averageplan<-c(fixhosp$Supply,rep ((NewCapacity/newH), newH ) )
  originalplan<-c(fixhosp$Supply,rep (0, newH ) )
  results<-data.frame(cbind(originalplan,averageplan,qp$solution,qp1$solution))
  results<-round(results,1)
  names(results)<-c('Original','Average','All','New')
  
  finalA<-data.frame(t(t(results)%*% t(Fij)))
  names(finalA)<-c('OrigAcc','AverAcc','AllAcc','NewAcc')
  
  
  print('Summary of the Standard deviation of differenct scenarioes')
  print(paste('Extant situation ',round(sd(finalA$OrigAcc),5)))
  print(paste('Assign average to new facilities ',round(sd(finalA$AverAcc),5))) 
  print(paste('MEA optimization for new facilities ',round(sd(finalA$NewAcc),5)) )
  print(paste('MEA optimization for All facilities ',round(sd(finalA$AllAcc),5))) 
  
  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Saving...")
  arc.progress_pos(90)
  
  Out_supply_data <- out_params[[1]]
  Out_Acc_data  <- out_params[[2]]
  
  arc.write(Out_supply_data, results,overwrite = TRUE)
  arc.write(Out_Acc_data, finalA,overwrite = TRUE)

  arc.progress_pos(100)
  }
