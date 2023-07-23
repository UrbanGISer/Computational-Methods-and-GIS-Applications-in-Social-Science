# CMGIS_Rtool: Computational Methods and GIS Applications in Social Sciences,Third Edition , CRC Press
# coding contributor: Lingbo Liu
# Center for Geographic Analysis, Harvard University
# School of Urban Design, Wuhan University

tool_exec <- function(in_params, out_params) {

  #####################################################################################################
  ### Check/Load Required Packages
  #####################################################################################################
  arc.progress_label("Loading packages...")
  arc.progress_pos(0)
  if (!requireNamespace("broom", quietly = TRUE)) {
    install.packages("broom", quiet = TRUE)
  }
  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    install.packages("minpack.lm", quiet = TRUE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    install.packages("tidyr", quiet = TRUE)
  }

  require(broom)
  require(tidyr)
  require(minpack.lm)

  arc.progress_label("Packages Load Successfully")
  arc.progress_pos(20)

  #####################################################################################################
  ### Define input/output parameters
  #####################################################################################################
  input_data <- in_params[[1]]
  dependent_variable <- in_params[[2]]
  independent_variables <- in_params[[3]]
  regressionmodel <- in_params[[4]] # Linear and NonLinear
  nonlinearfunction <- in_params[[5]]
  nonlinearstart <- in_params[[6]]
  weightmodel <- in_params[[7]]
  weightfield <- in_params[[8]]
  ExportorNot <- in_params[[9]]
  regressionsummary <- out_params[[1]]
  output_data <- out_params[[2]]
  FilterorNot <- in_params[[10]] # Regression with seleted records
  FilterCondition <- in_params[[11]]

  #####################################################################################################
  ### Load Data and Create Dataframe R Object
  #####################################################################################################

  arc.progress_label("Loading data...")
  arc.progress_pos(40)
  d <- arc.open(input_data)
  print("Read Data Successfully")
  df <- arc.select(d)
  # df <- arc.data2sf(dft)
  # df <- st_set_geometry(df, NULL)
  # Handle weighting field
  if (weightmodel == "On") {
    fields_list <- append(dependent_variable, c(independent_variables))
    fields_list <- append(fields_list, weightfield)
  } else {
    fields_list <- append(dependent_variable, c(independent_variables))
  }
  dfreg <- df[, names(df) %in% fields_list]

  # Handle Filter condition
  if (FilterorNot == "Yes") {
    df <- subset(df, eval(parse(text = FilterCondition)))
    dfreg <- df[, names(df) %in% fields_list]
    print(paste("Filter recoords with ", FilterCondition))
  }

  # Renamed data frame as Y,x1,x2,......, Wt
  xid <- 1:length(independent_variables)
  xnames <- paste("x", xid, sep = "")
  print("Orginal Variables:")
  print(names(dfreg))
  names(dfreg)[names(dfreg) == dependent_variable] <- "Y"
  names(dfreg)[names(dfreg) %in% c(independent_variables)] <- xnames
  if (weightmodel == "On") {
    names(dfreg)[names(dfreg) == weightfield] <- "wt"
  }
  print("Variables Renamed As :")
  print(names(dfreg))

  #####################################################################################################
  ### Define Linear Output Parameter
  #####################################################################################################

  LMresult <- function(regg) {
    # get coefficients
    tidyreg <- as.data.frame(tidy(regg))
    # round coefficients
    tidyreg[, -1] <- round(tidyreg[, -1], 4)
    # extract Name as row
    tidyname <- as.data.frame(t(as.data.frame(names(tidyreg))))
    # unify column names
    names(tidyreg) <- names(tidyname)
    # concatenate name and original data frame
    tidyreg <- rbind(tidyname, tidyreg)

    # get Statistics sumary
    glancereg <- as.data.frame(glance(regg))
    # extract important parameters
    glancereg <- glancereg[c(2, 5, 7:9)]
    # round coefficients
    glancereg <- round(glancereg, 4)
    # extract Name as row
    glancename <- as.data.frame(t(as.data.frame(names(glancereg))))
    # unify column names
    names(glancereg) <- names(glancename)
    # concatenate name and original data frame
    glancereg <- rbind(glancename, glancereg)
    # combine result
    all <- rbind(tidyreg, glancereg)
    # rename row names
    rownames(all) <- 1:nrow(all)
    return(all)
  }
  #####################################################################################################
  ### Define NoneLinear Output Parameter
  #####################################################################################################

  NLMresult <- function(regg) {
    # get coefficients
    tidyreg <- as.data.frame(tidy(regg))
    # round coefficients
    tidyreg[, -1] <- round(tidyreg[, -1], 4)
    # extract Name as row
    tidyname <- as.data.frame(t(as.data.frame(names(tidyreg))))
    # unify column names
    names(tidyreg) <- names(tidyname)
    # concatenate name and original data frame
    tidyreg <- rbind(tidyname, tidyreg)

    # Calculate Statistics summary
    RSS <- deviance(regg)
    R2 <- cor(dfreg$Y, fitted(regg))^2
    n <- nobs(regg)
    p <- length(coef(regg))
    adjR2 <- 1 - (1 - R2) * (n - 1) / (n - p - 1)
    pvalue <- coef(summary(regg))[1, 4]
    regvalue <- c(adjR2, pvalue, glance(regg)[4:6])
    regvalue <- lapply(regvalue, round, 4)
    regname <- c("adj.r.squared", "p.value", "logLik", "AIC", "BIC")
    glancereg <- data.frame(t(cbind(regname, regvalue)))
    # unify column names
    names(glancereg) <- names(tidyname)
    # combine result
    all <- rbind(tidyreg, glancereg)
    # rename row names
    rownames(all) <- 1:nrow(all)
    return(all)
  }
  #####################################################################################################
  ### Linear regression and Nonlinear Regression
  #####################################################################################################
  arc.progress_label("Run Regression")
  arc.progress_pos(60)

  RegTable <- data.frame(matrix(ncol = 5, nrow = 0))
  xdfbindname <- c("term", "estimate", "std.error", "statistic", "p.value")
  colnames(RegTable) <- xdfbindname

  if (regressionmodel == "Linear") {
    if (weightmodel == "On") {
      xnamesplus <- paste(xnames, collapse = "+")
      lmfunction <- paste("Y~", xnamesplus, sep = "")
      reglm <- lm(eval(parse(text = lmfunction)),
        data = dfreg,
        weights = wt
      )
    } else {
      f <- Y~ .
      reglm <- lm(f, data = dfreg)
    }
    print("Summary for Linear regression")
    print(summary(reglm))
    RegTable <- LMresult(reglm)
    df$Predict <- predict(reglm)
  } else if (regressionmodel == "Nonlinear") {
    startline <- paste("list(", nonlinearstart, ")", sep = "")
    functionline <- paste("Y~", nonlinearfunction, sep = "")
    if (weightmodel == "On") {
      reglm <- nlsLM(eval(parse(text = functionline)),
        data = dfreg,
        start = eval(parse(text = startline)),
        control = nls.lm.control(maxiter = 1000),
        weights = wt
      )
    } else {
      reglm <- nlsLM(eval(parse(text = functionline)),
        data = dfreg,
        start = eval(parse(text = startline)),
        control = nls.lm.control(maxiter = 1000)
      )
    }
    print("Summary for NoneLinear regression")
    print(summary(reglm))
    RegTable <- NLMresult(reglm)
    df$Predict <- predict(reglm)
  } else {
    for (i in xid) {
      lmfunction <- paste("Y~", xnames[i], sep = "")
      if (weightmodel == "On") {
        reglm <- lm(eval(parse(text = lmfunction)),
          data = dfreg,
          weights = wt
        )
      } else {
        reglm <- lm(eval(parse(text = lmfunction)),
          data = dfreg
        )
      }
      print("Summary for Linear regression")
      print(summary(reglm))
      RegTable1 <- LMresult(reglm)
      colnames(RegTable1) <- RegTable1[1, ]
      RegTable1 <- RegTable1[-1, ]
      row.names(RegTable1) <- NULL
      RegTable1[2, 1] <- independent_variables[i]
      RegTable <- rbind(RegTable, RegTable1)
    }
  }


  #####################################################################################################
  ### Write Output
  #####################################################################################################
  arc.progress_label("Print output...")
  arc.progress_pos(80)
  if (regressionmodel == "Linear") {
    colnames(RegTable) <- RegTable[1, ]
    RegTable <- RegTable[-1, ]
    row.names(RegTable) <- NULL
  }

  print("Statistical  Summary for  Regression")
  print(RegTable)
  print(paste("Summary saved as txt in : ", regressionsummary))
  if (regressionmodel == "Nonlinear") {
    RegTable <- as.data.frame(lapply(RegTable, as.character))
    # RegTable <- apply(RegTable, 1, as.character)
    write.csv(RegTable, regressionsummary)
  } else {
    arc.write(regressionsummary, data = RegTable, overwrite = TRUE)
  }
  if (ExportorNot == "Yes") {
    arc.write(output_data, data = df, overwrite = TRUE)
  }
  arc.progress_pos(100)
}
