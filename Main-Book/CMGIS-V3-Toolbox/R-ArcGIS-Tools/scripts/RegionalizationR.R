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
  
  if(!requireNamespace("rgeoda", quietly = TRUE))
    install.packages("rgeoda", quiet = TRUE) 
  if(!requireNamespace("sf", quietly = TRUE))
    install.packages("sf", quiet = TRUE) 

  require(rgeoda)
  require(sf)
  
  arc.progress_label("Packages Load Successfully")
  arc.progress_pos(20)
  
  ##################################################################################################### 
  ### Define input/output parameters
  #####################################################################################################
  input_data <- in_params[[1]]
  control_variables <- in_params[[2]]
  bound_variable <- in_params[[3]]
  bound_minValue<- in_params[[4]] 
  modeltype <- in_params[[5]] # REDcap SCHC AZP maxp
  redcap_param <- in_params[[6]] #fullorder-completelinkage
  schc_param <- in_params[[7]] # Complete
  number_cluster <- in_params[[8]]
  weight_model <- in_params[[9]] # queen or rook or knn
  knn_number <- in_params[[10]] #  knn
  output_data <- out_params[[1]]
  
  ##################################################################################################### 
  ### Load Data and Create Dataframe R Object 
  #####################################################################################################
  arc.progress_label("Loading data...")
  arc.progress_pos(40)
  d <- arc.open(input_data)
  print("Read Data Successfully")
  df <- arc.select(d)
  spdf <- arc.data2sf(df)
  
  # Prepare data
  control_list <- c(control_variables)
  # set control data
  clusterdata <-spdf[,names(spdf) %in% control_list ]
  # set bound variable
  bound_data <- spdf[bound_variable]

  ##################################################################################################### 
  ### Generate spatial weighs with specific model 
  #####################################################################################################
  arc.progress_label("Generating spatial weights...")
  arc.progress_pos(60)  
  # Handle weighting field
  if (weight_model=='Queen'){
    spatialweight <- queen_weights(spdf)
    } else if (weight_model=='Rook') {
      spatialweight <-  rook_weights(spdf)
    } else {
      spatialweight <-  knn_weights(spdf,knn_number)
    } 
  print(paste("Generate Spatial Weight with",weight_model,"Model",sep=""))
  
  ##################################################################################################### 
  ### Run Regionalization Model 
  #####################################################################################################
  arc.progress_label("Run Regionalization...")
  arc.progress_pos(80)
  # REDcap SCHC AZP maxp
  if (modeltype=='REDCAP'){
    reg_clusetr<- redcap(number_cluster, spatialweight, clusterdata, redcap_param,bound_variable =bound_data ,min_bound = bound_minValue)
  } else if (modeltype=='SCHC'){
    reg_clusetr<- schc(number_cluster, spatialweight, clusterdata, schc_param, bound_variable =bound_data ,min_bound = bound_minValue)
  } else if (modeltype=='AZP_greedy'){
    reg_clusetr<- azp_greedy(number_cluster, spatialweight, clusterdata,bound_variable =bound_data ,min_bound = bound_minValue)
  } else if (modeltype=='SKATER'){
    reg_clusetr<- skater(number_cluster, spatialweight, clusterdata,bound_variable =bound_data ,min_bound = bound_minValue)    
  } else {
    reg_clusetr<- maxp_greedy(spatialweight, clusterdata,bound_variable =bound_data ,min_bound = bound_minValue)
  }
  print("Regionalization Finished")
  print (reg_clusetr)  

  
  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Saving...")
  arc.progress_pos(90)
  
  df$cluster<-reg_clusetr$Clusters
  arc.write(output_data, df,overwrite = TRUE, shape_info = arc.shapeinfo(d))
  #print("Regresion summary table saved")
  arc.progress_pos(100)
  }
