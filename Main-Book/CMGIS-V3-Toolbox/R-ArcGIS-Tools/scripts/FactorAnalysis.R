# CMGIS_Rtool: Computational Methods and GIS Applications in Social Sciences,Third Edition , CRC Press
# Principal Component Analysis and Factor Analysis
tool_exec <- function(in_params, out_params) {

  #####################################################################################################
  ### Check/Load Required Packages
  #####################################################################################################
  arc.progress_label("Loading packages...")
  arc.progress_pos(20)

  if (!requireNamespace("psych", quietly = TRUE)) {
    install.packages("psych", quiet = TRUE)
  }
  if (!requireNamespace("GPArotation", quietly = TRUE)) {
    install.packages("GPArotation", quiet = TRUE)
  }
  require(psych)
  require(GPArotation)
  #####################################################################################################
  ### Define input/output parameters
  #####################################################################################################
  input_data <- in_params[[1]]
  excluded_variables <- in_params[[2]]
  n_factor <- (in_params[[3]])
  modelchoose <- in_params[[4]]
  FA_rotate <- in_params[[5]]
  FA_fm <- in_params[[6]]
  PCA_rotate <- in_params[[7]]
  plothow <- in_params[[8]]
  sortloading <- in_params[[9]]
  output_prediction_data <- out_params[[1]]

  #####################################################################################################
  ### Load Data and Create Dataframe R Object
  #####################################################################################################
  arc.progress_label("Loading data...")
  arc.progress_pos(40)

  d <- arc.open(input_data)
  exclued_list <- c(excluded_variables)
  fields_list <- names(d@fields)
  df0 <- arc.select(d, fields = fields_list)
  df_keep <- df0[, !names(df0) %in% exclued_list]
  df_exclude <- df0[, names(df0) %in% exclued_list]
  df <- df_keep
  #####################################################################################################
  ### Factor Analysis
  #####################################################################################################
  arc.progress_label("Creating training and testing datasets...")
  arc.progress_pos(60)

  if (modelchoose == "Factor Analysis") {
    result_pca_fa <- fa(df, nfactors = n_factor, rotate = FA_rotate, fm = FA_fm)
  } else {
    result_pca_fa <- principal(df, nfactors = n_factor, rotate = PCA_rotate, score = TRUE)
  }
  print(result_pca_fa, sort = sortloading, digits = 4)
  cat("\n ---------------------------\n Eigen Values of all components\n")
  # get the eigenvalue and Proportion
  r1 <- result_pca_fa$values
  r2 <- r1 / sum(r1)
  r3 <- cumsum(r2)
  dfr13 <- data.frame(r1, r2, r3)
  names(dfr13) <- c("Eigenvalue", "Proportion", "Cumulative")
  print(round(dfr13, digits = 4))

  result <- as.data.frame(result_pca_fa$scores)
  final_result <- cbind(df_exclude, result)
  names(final_result)[names(final_result) == "df_exclude"] <- exclued_list
  # Scree plots with parallel analysis
  if (plothow == "None") {
    print("Igenor all plots ")
  } else if (plothow == "Scree plot") {
    print("Scree plots with parallel analysis")
    fa.parallel(df, n.obs = 221, fa = "both", n.iter = 100, main = "Scree plots with parallel analysis")
  } else {
    print("Factor loading matrices")
    fa.diagram(result_pca_fa)
  }


  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Writing output...")
  arc.progress_pos(80)

  arc.write(output_prediction_data, final_result, overwrite = TRUE)

  arc.progress_pos(100)
}
