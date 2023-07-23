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
  
  if(!requireNamespace("lpSolve", quietly = TRUE))
    install.packages("lpSolve", quiet = TRUE) 
  if(!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", quiet = TRUE) 

  require(lpSolve)
  require(dplyr)
  
  arc.progress_label("Packages Load Successfully")
  arc.progress_pos(10)
  
  ##################################################################################################### 
  ### Define input/output parameters
  #####################################################################################################
  input_data <- in_params[[1]]
  OriginID<- in_params[[2]]
  Work <- in_params[[3]]
  DestinationID<- in_params[[4]] 
  Emp <- in_params[[5]] # REDcap SCHC AZP maxp
  Netwtime <- in_params[[6]] #fullorder-completelinkage
  output_data <- out_params[[1]]
  
  ##################################################################################################### 
  ### Load Data and Create Dataframe R Object 
  #####################################################################################################
  arc.progress_label("Loading data...")
  arc.progress_pos(30)
  d <- arc.open(input_data)
  print("Read Data Successfully")
  df <- arc.select(d)

  # Prepare data
  fields_list=c(OriginID,Work,DestinationID,Emp,Netwtime)
  # Choose selected data
  data_od<- subset(df, select = fields_list)
  # Standardize column names
  names(data_od)<-c('OID','WORK','DID','EMP','NetwTime')
  #Sort data first based on origin TAZ code, second based on destination TAZ code
  data_od<-data_od[order(data_od[,"OID"],data_od[,"DID"]),]

  ##################################################################################################### 
  ### Prepare Matrix data for subsequent analysis
  #####################################################################################################
  arc.progress_label("Building Matrix...")
  arc.progress_pos(50) 
  
  #Then convert the OD time to a vector
  vec_time<-data_od[,"NetwTime"]
  vec_time<-t(vec_time)
  print( paste('The average commute time  is ',round(mean(vec_time),3)))
  
  
  #Read resident works constraints into a data frame
  data = data_od %>% group_by(OID) %>%
    summarise(WORK = first(WORK),
              .groups = 'drop')
  #Sort data based on origin TAZ code
  data<-data[order(data[,"OID"]),]
  #Convert the number of resident workers to vectors
  vec_res<-data[,"WORK"]
  vec_res<-t(vec_res)
  print( paste('The population of total commuters is ',sum(vec_res)))
  
  #Read job constraints  into a data frame
  data = data_od %>% group_by(DID) %>%
    summarise(EMP = first(EMP),
              .groups = 'drop')
  #Sort data based on destination TAZ code
  data<-data[order(data[,"DID"]),]
  #Convert the number of jobs to vectors
  vec_emp<-data[,"EMP"]
  vec_emp<-t(vec_emp)
  
  #Calculate # of REZ and EMP
  nREZ=length(vec_res)
  nEMP=length(vec_emp)
  
  #Convert the above vector to a matrix for later use
  costs<-matrix(vec_time,nrow=nREZ,ncol=nEMP,byrow=T)
  
  #Set up constraint signs and right-hand sides
  row.signs<-rep("=",nREZ)
  row.rhs<-vec_res
  col.signs<-rep("<=",nEMP)
  col.rhs<-vec_emp
  
  ##################################################################################################### 
  ### Run lpSolve package
  #####################################################################################################
  arc.progress_label("Apply lpSolve...")
  arc.progress_pos(60)
  #Run to measure wasteful commuting
  lpresult<-lp.transport(costs,"min",row.signs,row.rhs,col.signs,col.rhs)
  print( paste('Value of objective function at optimum is ',lpresult$objval))
  print( paste('The average optimum commute time  is ',round((lpresult$objval/sum(vec_res)),3)))
  wastc<-round((mean(vec_time)-lpresult$objval/sum(vec_res))/mean(vec_time)*100,3)
  print( paste('wasteful commuting is ',wastc,'%'))
  
  # Result to dataframe
  result<-data.frame(data_od,as.vector(t(lpresult$solution)),row.names=NULL)
  colnames(result)<-c('OID','WORK','DID','EMP','NetwTime','Flow')
  
  
  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Saving...")
  arc.progress_pos(90)
  
  arc.write(output_data, result,overwrite = TRUE)

  #print("Regresion summary table saved")
  arc.progress_pos(100)
  }
