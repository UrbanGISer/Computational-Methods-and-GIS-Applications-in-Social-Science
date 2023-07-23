# CMGIS_Rtool: Computational Methods and GIS Applications in Social Sciences,Third Edition , CRC Press
# Principal Component Analysis and Factor Analysis
tool_exec<- function(in_params, out_params){
  
  #####################################################################################################  
  ### Check/Load Required Packages  
  #####################################################################################################   
  arc.progress_label("Loading packages...")
  arc.progress_pos(80)
  
  if(!requireNamespace("spatialreg", quietly = TRUE))
    install.packages("spatialreg", quiet = TRUE)
  if(!requireNamespace("spdep", quietly = TRUE))
    install.packages("spdep", quiet = TRUE)
  if(!requireNamespace("broom", quietly = TRUE))
    install.packages("broom", quiet = TRUE) 
  if(!requireNamespace("sp", quietly = TRUE))
    install.packages("sp", quiet = TRUE) 
  if(!requireNamespace("tidyr", quietly = TRUE))
    install.packages("tidyr", quiet = TRUE) 
  require(spatialreg)
  require(broom)
  require(sp)
  require(tidyr)
  require(spdep)
  ##################################################################################################### 
  ### Define input/output parameters
  #####################################################################################################
  input_data <- in_params[[1]]
  dependent_variable <- in_params[[2]]
  independent_variables <- in_params[[3]]
  SWMmodel<- in_params[[4]] #FastQuenn,Custummized_SWM_Table
  inputSWMtable <- in_params[[5]]
  output_prediction_data <- out_params[[1]]

  ##################################################################################################### 
  ### Load Data and Create Dataframe R Object 
  #####################################################################################################
  arc.progress_label("Loading data...")
  arc.progress_pos(80)
  
  d <- arc.open(input_data)
  fields_list <- append(c(dependent_variable), independent_variables)
  df <- arc.select(d)
  dfreg <- arc.select(d, fields = fields_list)
  spdf <- arc.data2sp(df)
  spdfreg <- arc.data2sp(dfreg)
  yname<-dependent_variable
  #####################################################################################################
  ### Build matrix
  #####################################################################################################
  arc.progress_label("Building Spatial Weight Matrix...")
  arc.progress_pos(80)
  
  if (SWMmodel=="Queen") {
    nbf <- poly2nb(spdf, queen=T)
    listweight<- nb2listw(nbf, style="W", zero.policy=TRUE)
  } else {
    SWM <- arc.open(inputSWMtable)
    SWMtable <- arc.select(SWM)
    SWMtable <-SWMtable[-1]
    names(SWMtable)<-c("O","D","w")
    SWMtableF <-spread(SWMtable, key=D, value=w)
    row.names(SWMtableF)<-SWMtableF[,1]
    SWMtableF <-SWMtableF[-1]
    SWMtableF[is.na(SWMtableF)] = 0
    listweight<- mat2listw(as.matrix(SWMtableF), style="W")
  }
  
  #####################################################################################################
  ### OLS regression and Spatial Autocorrelation Test
  #####################################################################################################
  arc.progress_label("OLS regression and Spatial Autocorrelation Test...")
  arc.progress_pos(80)
  
  names(spdfreg)[names(spdfreg) == yname] <- "Y"
  f<-Y ~ .
  print("#############################################")
  print("OLS regression")
  reglm=lm(f,data=spdfreg)
  print(summary(reglm))
  
  # call spatial test functions
  print("#############################################")
  print("Moran's I test")
  print(lm.morantest(reglm,listweight))
  print("#############################################")
  print("LM test")
  print(lm.LMtests(reglm,listweight,test=c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")))
  
  #####################################################################################################
  ### RUN Spatial Regression Model
  #####################################################################################################
  arc.progress_label("Running Spatial Regression Model...")
  arc.progress_pos(80)
  
  #SAR Spatial Lag (Autoregressive) Model y=pWy+XB+e 
  SLMreg=lagsarlm(f,data= spdfreg, listweight)
  print("#############################################")
  print("Spatial Lag Model")  
  print(summary(SLMreg,Nagelkerke=T))
  # call more functions
  #SEM Spatial Error Model  y=XB+u,   u=LWu+e
  SEMreg=errorsarlm(f,data=spdfreg, listweight)
  print("#############################################")
  print("Spatial Error Model")  
  print(summary(SEMreg,Nagelkerke=T))
  
  #####################################################################################################
  ### Collect Result as Table
  #####################################################################################################
  arc.progress_label("Preparing Regression Result...")
  arc.progress_pos(80)
  
  Regresult <- function(regg,modelname) {
    tidyreg <- as.data.frame(tidy(regg))
    tidyreg$model<-modelname
    tidyname<-as.data.frame(t(as.data.frame(names(tidyreg))))
    rownames(tidyname)<-modelname
    names(tidyreg)<-names(tidyname)
    rownames(tidyreg)<-1:nrow(tidyreg)
    tidyreg<-rbind(tidyname,tidyreg)
    #glance
    glancereg<- as.data.frame(glance(regg))
    glancename<-as.data.frame(t(as.data.frame(names(glancereg))))
    names(glancereg)<-names(glancename)
    rownames(glancename)<-"Fitness"
    rownames(glancereg)<-"Value"
    glancereg<-rbind(glancename,glancereg)
    #Combine
    all<-rbind(tidyreg,glancereg)
    return(all)
  }
  
  SLMtable<-Regresult(SLMreg,"SLM")
  SEMtable<-Regresult(SEMreg,"SEM")
  Alltable<-rbind(SLMtable,SEMtable)
  names(Alltable)<-SLMtable[1,]
  
  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Writing output...")
  arc.progress_pos(20)
  
  arc.write(output_prediction_data, Alltable,overwrite = TRUE)
  arc.progress_pos(100)
}
