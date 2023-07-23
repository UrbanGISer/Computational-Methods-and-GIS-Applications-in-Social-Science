# CMGIS_Rtool: Computational Methods and GIS Applications
# in Social Sciences,Third Edition , CRC Press
# Code contributor : Lingbo Liu
# Center for Geographic Analysis, Harvard University
# School of Urban Design, Wuhan University

tool_exec <- function(in_params, out_params) {

    ####################################################
    ### Check/Load Required Packages
    ####################################################
    arc.progress_label("Loading packages...")
    arc.progress_pos(0)

    if (!requireNamespace("stats", quietly = TRUE)) {
        install.packages("stats", quiet = TRUE)
    }
    suppressPackageStartupMessages({
        require(stats)
        require(factoextra)
    })

    arc.progress_label("Packages Load Successfully")
    arc.progress_pos(10)

    ##################################################
    ### Define input/output parameters
    ##################################################
    # Input data
    input_data <- in_params[[1]]
    cluster_variable <- in_params[[2]]
    ncluster <- in_params[[3]]
    # Simulation times
    clmodel <- in_params[[4]]
    clmethod <- in_params[[5]]
    # Output Estimate
    plotornot <- in_params[[6]]
    plotoptimal <- in_params[[7]]
    output_data <- out_params[[1]]

    ##################################################
    ### Load Data and Create Dataframe R Object
    ##################################################
    arc.progress_label("Loading data...")
    arc.progress_pos(20)

    d <- arc.open(input_data)
    print("Read Data Successfully")
    df <- arc.select(d)
    fields_list <- c(cluster_variable)
    dfreg <- df[, names(df) %in% fields_list]
    print("Data chose Successfully")
    # total within sum of square
    TSS <- function(x, g) {
        sum(aggregate(x, by = list(g), function(x) {
            sum(scale(x,
                scale = FALSE
            )^2)
        })[, -1])
    }

    if (clmodel == "k-means") {
        km <- kmeans(dfreg, centers = ncluster, nstart = 5)
        df$clusterID <- factor(km$cluster)
        if (plotoptimal == "Yes") {
            wss <- (nrow(dfreg) - 1) * sum(apply(dfreg, 2, var))
            for (i in 2:20) wss[i] <- sum(kmeans(dfreg, centers = i)$withinss)
            clusternumber <- 1:20
            WithinSumSquare <- wss
            plot(1:20, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
            print(data.frame(clusternumber, WithinSumSquare))
        }
    } else {
        dfdist <- dist(dfreg)
        hc <- hclust(dfdist, method = clmethod)
        hcclusters <- cutree(hc, k = ncluster)
        df$clusterID <- factor(hcclusters)
        if (plotornot == "Yes") {
            plot(hc)
        }
        if (plotoptimal == "Yes") {
            x.grps <- cutree(hc, 1:20)
            TSS.all <- apply(x.grps, 2, function(g) TSS(dfreg, g))
            plot(1:20, TSS.all, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
            clusternumber <- 1:20
            WithinSumSquare <- TSS.all
            print(data.frame(clusternumber, WithinSumSquare))
        }
    }

    ##################################################
    ### Export data
    ##################################################
    arc.progress_label("Saving...")
    arc.progress_pos(90)

    arc.write(output_data, data = df, overwrite = TRUE)

    arc.progress_pos(100)
}
