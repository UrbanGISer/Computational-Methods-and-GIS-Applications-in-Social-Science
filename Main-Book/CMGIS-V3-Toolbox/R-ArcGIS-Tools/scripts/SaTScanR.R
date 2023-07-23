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

    if (!requireNamespace("sf", quietly = TRUE)) {
        install.packages("sf", quiet = TRUE)
    }
    if (!requireNamespace("spatstat", quietly = TRUE)) {
        install.packages("spatstat", quiet = TRUE)
    }
    if (!requireNamespace("smacpod", quietly = TRUE)) {
        install.packages("smacpod", quiet = TRUE)
    }
    if (!requireNamespace("maptools", quietly = TRUE)) {
        install.packages("maptools", quiet = TRUE)
    }
    suppressPackageStartupMessages({
        require(sf, quietly = T)
        require(spatstat, quietly = T)
        require(smacpod, quietly = T)
        require(maptools, quietly = T)
    })
    arc.progress_label("Packages Load Successfully")
    arc.progress_pos(10)

    ##################################################
    ### Define input/output parameters
    ##################################################
    # Input data
    control_feature <- in_params[[1]]
    case_feature <- in_params[[2]]
    boundary_feature <- in_params[[3]]
    # Simulation times
    nsim <- in_params[[4]]
    # Output Estimate
    case_result <- out_params[[1]] # estimate gtif
    circle_result <- out_params[[2]] # test data


    ##################################################
    ### Load Data and Create Dataframe R Object
    ##################################################
    arc.progress_label("Loading data...")
    arc.progress_pos(20)

    arc_to_sf <- function(inputfeautre) {
        d <- arc.open(inputfeautre)
        df <- arc.select(d)
        spdf <- arc.data2sf(df)
        return(spdf)
    }
    arc_to_sfpt <- function(inputfeautre, category) {
        spdf <- arc_to_sf(inputfeautre)
        spt <- data.frame(st_coordinates(spdf))
        spt$marks <- category
        return(spt)
    }


    # Get control and case point
    ctrl_sp <- arc_to_sfpt(control_feature, "control")
    case_sp <- arc_to_sfpt(case_feature, "case")
    # For output
    casept <- arc_to_sf(case_feature)
    # Combination data
    df <- rbind(ctrl_sp, case_sp)
    rownames(df) <- NULL

    # Get Boundary for window
    # must import maptools, otherwise owin fail
    boundary <- arc_to_sf(boundary_feature)

    boundary_owin <- as.owin(as_Spatial(boundary))

    ##################################################
    ### Build PPP data
    ##################################################
    arc.progress_label("Build PPP ...")
    arc.progress_pos(30)

    px <- ppp(df$X, df$Y,
        window = boundary_owin,
        marks = df$marks
    )

    ##################################################
    ### Spatial Scan Test
    ##################################################
    arc.progress_label("Apply Spatial Scan Test...")
    arc.progress_pos(50)

    scan <- spscan.test(px, nsim = nsim, case = "case")

    # plot scan test results
    # plot(scan, chars = c(4, 20), main = "most likely cluster for grave data",
    #     border = "red")

    # Extact result
    clu_circles <- st_sf(st_sfc())
    clu_id <- data.frame()
    for (i in (1:length(scan$clusters))) {
        # Generate Circle
        lon <- scan$clusters[[i]]$coords[1, 1]
        lat <- scan$clusters[[i]]$coords[1, 2]
        r <- scan$clusters[[i]]$r
        tower <- st_point(x = c(lon, lat), dim = "XY")
        clu_center <- st_point(x = c(lon, lat), dim = "XY")
        dat_circles <- st_buffer(clu_center, dist = r)
        datatemp <- data.frame(clusterID = i)
        datatemp$geom <- st_as_text(st_geometry(dat_circles))
        datatemp <- st_as_sf(datatemp, wkt = "geom")
        datatemp$loglikrat <- scan$clusters[[i]]$loglikrat
        datatemp$pvalue <- scan$clusters[[i]]$pvalue
        datatemp$centerx <- scan$clusters[[i]]$coords[1, 1]
        datatemp$centery <- scan$clusters[[i]]$coords[1, 2]
        datatemp$r <- scan$clusters[[i]]$r
        clu_circles <- rbind(clu_circles, datatemp)
        # Get ID
        datatemp <- data.frame(pointID = scan$clusters[[i]]$locids)
        datatemp$clusterID <- i
        clu_id <- rbind(clu_id, datatemp)
    }

    # Set CRS information
    st_crs(clu_circles) <- st_crs(casept)$input
    clu_case <- clu_id[clu_id$pointID > nrow(ctrl_sp), ]
    clu_case$pointID <- clu_case$pointID - nrow(ctrl_sp)

    casept$pointID <- 1:nrow(casept)
    casept <- merge(casept, clu_case, by = "pointID", all.x = TRUE)

    # kd <-kdest(px, case = case = "case")
    # plot(kd, cbind(iso, theo) ~ r, legend = FALSE, main = "")
    # kdenv <- kdest(px case = "affected", nsim = 49, level = 0.95)
    # print(kdenv)
    # plot(kdenv, ylab = "difference in K functions")
    # legend("topleft", legend = c("obs", "avg", "max/min env", "95% env"),
    #    lty = c(1, 2, 1, 2), col = c("black", "red", "darkgrey", "lightgrey"),
    #    lwd = c(1, 1, 10, 10))


    ##################################################
    ### Export data
    ##################################################
    arc.progress_label("Saving...")
    arc.progress_pos(90)

    arc.write(case_result, data = casept, overwrite = TRUE)
    arc.write(circle_result, data = clu_circles, overwrite = TRUE)

    arc.progress_pos(100)
}
