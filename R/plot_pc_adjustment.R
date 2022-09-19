#' @title Plot MCN Performance Against Number of PCs Included.
#' @description Plot the correlation between estimated MCN and external phenotype (or MCN PGS if no phenotype is provided) against number of PCs used for adjustment.
#' @param output_dir Path to the output directory of `plot_pc_adjustment` function with "tempdir" directory in it.
#' @param max_pc_num (Optional) A number indicating the maximum principal components to be used for plotting. Either a fraction between 0 and 1 indicating proportion of the sample size or an positive integer indicating count of PCs are accepted. Default is `0.01`.
#' @param correlated_pheno (Optional) A data frame with columns "Sample_Name" and "Phenotype". The phenotype should be associated with MCN for checking the minimum number of PCs needed to capture the confounding factors. Default is `NULL` which means the MCN PGS will be used for PC checking. Note that the association between MCN PGS and MCN estimates may not strong enough in small populations (e.g. less than a few hundred of samples) so an associated phenotype is highly recommended.
#' @return A plot showing correlation between estimated MCN and external phenotype (or MCN PGS if no phenotype is provided) against number of PCs used for adjustment.
#' @export
#' @import dplyr
#' @import doParallel
#' @import foreach
#' @import logging
#' @importFrom parallel detectCores
#' @importFrom data.table fread setDTthreads
#' @importFrom HardyWeinberg HWExactStats
#' @importFrom stats cor lm prcomp
plot_pc_adjustment <- function(output_dir, max_pc_num = 0.01, correlated_pheno = NULL) {
    mcn_dir <- file.path(output_dir, "tempdir", "3_R_MCN")
    auto_path <- file.path(mcn_dir, "3_autosomal_LRR.rds")
    mito_path <- file.path(mcn_dir, "3_mitochondrial_LRR.rds")
    autopc_path <- file.path(mcn_dir, "4_autosomal_PCA.rds")
    mcnpgs_path <- file.path(mcn_dir, "7_pc_orient", "2_mcn_pgs.rds")
    if (!dir.exists(output_dir)) stop("Directory ", output_dir, " does not exist.")
    if (file.exists(file.path(output_dir, "check_pc_number.log"))) invisible(file.remove(file.path(output_dir, "check_pc_number.log")))
    logReset()
    addHandler(writeToConsole)
    addHandler(writeToFile, file = file.path(output_dir, "check_pc_number.log"))
    logger <- getLogger()
    logger$setLevel('DEBUG')
    setLevel('INFO', logger$handlers$writeToConsole)
    setLevel('DEBUG', logger$handlers$writeToFile)
    loginfo("Starts PC number check.")
    loginfo(paste0("Log file created at ", file.path(output_dir, "check_pc_number.log"), "."))
    loginfo("Checking arguments...")
    if (!file.exists(auto_path)) {
        logerror(paste0(auto_path, " does not exist. Please run estimate_mcn_from_cel() with \"keep_tempfile = TRUE\" to generate the proper intermediate results for PC number check."))
        stop(auto_path, "does not exist.")
    }
    if (!file.exists(mito_path)) {
        logerror(paste0(mito_path, " does not exist. Please run estimate_mcn_from_cel() with \"keep_tempfile = TRUE\" to generate the proper intermediate results for PC number check."))
        stop(mito_path, "does not exist.")
    }
    if (max_pc_num <= 0 | (max_pc_num >= 1 & as.integer(max_pc_num) != max_pc_num)) {
        logerror("\"max_pc_num\" should be a fraction between 0 and 1 (proportion of post-QC sample size) or an positive integer (count of PCs to be used).")
        stop("Invalid \"max_pc_num\" argument.")
    }
    mito_lrr <- readRDS(mito_path)
    if (is.null(correlated_pheno)) {
        if (ncol(mito_lrr) >= 1000) {
            loginfo("No correlated phenotype provided. Will use MCN PRS for PC number check.")
        } else {
            logwarn("No correlated phenotype provided. Will use MCN PRS for PC number check. Note that this may not powerful enough in small sample size (e.g. <1000).")
        }
        if (!file.exists(mcnpgs_path)) {
            logerror(paste0(mcnpgs_path, " does not exist. Please run estimate_mcn_from_cel() with \"keep_tempfile = TRUE\" and \"correlated_pheno = NULL\" to generate the proper intermediate results for PC number check."))
            stop(mcnpgs_path, "does not exist.")
        }
        correlated_pheno <- readRDS(mcnpgs_path) %>%
            rename(Phenotype = PRS)
        pheno_name <- "MCN_PGS"
    } else {
        pheno_name <- "Phenotype"
        if (!all(c("Sample_Name", "Phenotype") %in% colnames(correlated_pheno))) {
            logerror("\"correlated_pheno\" should be either NULL or a data frame with \"Sample_Name\" and \"Phenotype\" columns.")
            stop("Invalid \"correlated_pheno\" argument.")
        }
    }
    loginfo("All preparation complete.")
    if (max_pc_num < 1) {
        numpc <- round(max_pc_num * ncol(mito_lrr))
        if (numpc == 0) numpc <- 1
    } else {
        numpc <- max_pc_num
    }
    if (max_pc_num < 1 & numpc < 15) {
        logwarn(paste0("Only ", numpc, "PC(s) were used (", ncol(mito_lrr), "*", max_pc_num, "). Consider set a higher \"max_pc_num\" to better evaluate PC numbers to be used."))
    }
    if (numpc > ncol(mito_lrr)) {
        logwarn(paste0("Number of PCs used (", numpc, ") is larger than the post QC sample size (", ncol(mito_lrr), "). Will use ", ncol(mito_lrr), " PCs."))
        numpc <- ncol(mito_lrr)
    }
    flag_calcPC <- TRUE
    if (file.exists(autopc_path)) {
        auto_pca <- readRDS(autopc_path)
        if (ncol(auto_pca) >= numpc) {
            loginfo("Autosome PCA file already exist. Skipping PCA...")
            auto_pca <- auto_pca[, 1:numpc]
            flag_calcPC <- FALSE
        } else {
            loginfo(paste0("Existing autosome PCA file doesn't have enough PCs (", ncol(auto_pca), "<", numpc, "). Recalculating ", numpc, " PC(s)..."))
        }
    } else {
        loginfo(paste0("Calculating ", numpc, " PC(s)..."))
    }
    if (flag_calcPC) {
        auto_lrr <- readRDS(auto_path)
        auto_pca <- prcomp(t(auto_lrr), rank. = numpc, scale. = TRUE)$x
        logdebug("Saving autosomal PC(s)...")
        saveRDS(auto_pca, autopc_path)
    }
    loginfo("Calculating MCN estimates adjusted by different number of autosomal PCs...")
    registerDoParallel(detectCores())
    temp_res <- foreach(i = 1:numpc, .combine = c) %dopar% {
        mito_lrr_adjusted <- apply(mito_lrr, 1, function(X) {
            lm(X ~ auto_pca[, 1:i])$residuals
        })
        temp_merge <- prcomp(mito_lrr_adjusted, rank. = 1, scale. = TRUE)$x %>%
            as.data.frame() %>%
            rename(MCN = PC1) %>%
            mutate(Sample_Name = rownames(mito_lrr_adjusted)) %>%
            select(Sample_Name, MCN) %>%
            merge(correlated_pheno)
        cor(temp_merge$Phenotype, temp_merge$MCN, method = "spearman")
    }
    stopImplicitCluster()
    plot(abs(temp_res), xlab = "Number of PC(s) included", ylab = pheno_name)
}