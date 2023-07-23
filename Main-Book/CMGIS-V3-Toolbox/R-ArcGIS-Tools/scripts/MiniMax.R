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
  
  if(!requireNamespace("Rglpk", quietly = TRUE))
    install.packages("Rglpk", quiet = TRUE) 
  if(!requireNamespace("slam", quietly = TRUE))
    install.packages("slam", quiet = TRUE) 
  require(Rglpk)
  require(lpSolve)
  arc.progress_label("Packages Load Successfully")
  arc.progress_pos(10)
  
  ##################################################################################################### 
  ### Define input/output parameters
  #####################################################################################################
  input_data <- in_params[[1]]
  OriginID<- in_params[[2]]
  CandidatesList <- in_params[[3]]    # Distance matrix of Candidate Locations 
  Requiredchoice<-in_params[[4]]  # Yes or No
  RequiredColumn<- in_params[[5]] # Set the Nearest table as Last column
  k <- in_params[[6]] 
  output_data <- out_params[[1]]
  
  ##################################################################################################### 
  ### Load Data and Create Dataframe R Object 
  #####################################################################################################
  arc.progress_label("Loading data...")
  arc.progress_pos(30)
  d <- arc.open(input_data)
  all_list <-names(d@fields)
  df <- arc.select(d, fields = all_list )
  data_oid<- subset(df, select = c(OriginID))
  data_od<-df[,names(df) %in% CandidatesList ]
  filelist<-unlist(CandidatesList)
  if (Requiredchoice == "Yes"){
    data_rc<-subset(df, select = c(RequiredColumn))
    data_od<-cbind(data_od,data_rc)
    filelist<-c(filelist,RequiredColumn)
  } 
  
  ##################################################################################################### 
  ### Prepare Matrix data for subsequent analysis
  #####################################################################################################
  arc.progress_label("Building Matrix...")
  arc.progress_pos(50) 
  
  
  m<-ncol(data_od) # supply facilities
  n<-nrow(data_od) # demand locations
  
  # 01 a numeric vector representing the objective coefficients for MiniMax function min(z)
  coef <- c(rep (0 , m*(n+1) ) , 1)

  # 02 a numeric vector representing the right hand side of the constraints
  rhs <- c(rep (1 , n ) , rep (0 , n*(m+1)),k )
  
  # 03 binary or continuous value 
  var_type <- c(rep("B", (n+1)*m ) , "C")
  
  # 04 a numeric vector or a (sparse) matrix of constraint coefficients
  Amatrix <- matrix (0 ,n*(m+2)+1, length(coef) )

  # 04-1 SIGMA( Xij) =1  for i =1 ,, n
  for ( i in 1: n ) {
    for ( j in 1: m ) {
      Amatrix [i , m*(i -1) + j ] <-1
    }
  }
  # 04-2 SIGMA( Cij xij )<-z  for i =1 ,, n (assign Cijxij ) 
  for ( i in 1: n ) {
    for ( j in 1: m ) {
      Amatrix [(n+i) , m*(i -1) +j ] <- data_od[i , j ]
    }
  }

  # 04-2  SIGMA( Cij xij )<-z  for i =1 ,, n (assign -1 to z)
  Amatrix[(n+1):(2*n),length(coef)]<--1
  
  # 04-3 xij< yi (assign xij ) 
  for ( i in 1:(n*m)) {
    Amatrix [ (2*n+i) , i ] <- 1
  }
  # 04-3 xij< yi; (assign -yi) 
  for ( i in 1: n ) {
    for ( j in 1: m ) {
      Amatrix [(2*n+(i-1)*m+j) , (n*m +j) ] <- -1
    }
  } 

  # 04-4 SIGMA(Y) =k
  Amatrix [dim(Amatrix)[1],(n*m +1): (n*m +m) ] <- 1
  
  # 05 a character vector with the directions of the constraints
  signs <- c( rep("==", n ) , rep("<=", n*(m+1) ),"==" )
  

  # Prepare data
  if(Requiredchoice=='Yes' ){
    rhs <- c(rhs,1 )
    Amatrixlast<-c(rep (0 , (m*(n+1)-1) ) , 1,0)
    Amatrix  <-rbind(Amatrix ,Amatrixlast) 
    signs <- c( signs,  "==")
  } 

  ##################################################################################################### 
  ### Run Rglpk package
  #####################################################################################################
  arc.progress_label("Apply Rglpk...")
  arc.progress_pos(60)
  # lp package
  #lpresult <- lp( "min",coef ,Amatrix , signs , 
  #                rhs, binary.vec = 1:(length(coef)-1), compute.sens=TRUE)
  lpresult <- Rglpk_solve_LP ( obj =coef , mat=Amatrix , 
                               dir  = signs , types = var_type , 
                              rhs = rhs , max= FALSE )

  print( 'columns input  as ')
  print( filelist)
  print( 'columns status as ')
  print( lpresult$ solution[(n*m+1):(n*m+m)])
  print( paste('The Optimal minimum Maximun value is ',round(lpresult$ solution[length(coef)],3)))
  

  # Result to data frame
  result<-data.frame(matrix (lpresult$ solution[1:(n*m)]  , n , m , byrow = TRUE ))
  colnames(result)<-filelist
  resultfinal<-cbind(data_oid,result)
  
  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Saving...")
  arc.progress_pos(90)
  
  arc.write(output_data, resultfinal,overwrite = TRUE)

  arc.progress_pos(100)
  }
