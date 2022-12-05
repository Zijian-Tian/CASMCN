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
#' @import bigstatsr
#' @importFrom parallel detectCores
#' @importFrom data.table fread setDTthreads
#' @importFrom HardyWeinberg HWExactStats
#' @importFrom stats cor lm predict
#' @importFrom tidyr drop_na
plot_pc_adjustment <- function(output_dir, max_pc_num = 0.01, correlated_pheno = NULL) {
    mcn_dir <- file.path(output_dir, "tempdir", "3_R_MCN")
    auto_path <- file.path(mcn_dir, "3_autosomal_LRR.rds")
    mito_path <- file.path(mcn_dir, "3_mitochondrial_LRR.rds")
    info_path <- file.path(mcn_dir, "3_sample_QC.rds")
    autopc_path <- file.path(mcn_dir, "4_autosomal_PCA.rds")
    mcnpgs_path <- file.path(mcn_dir, "7_pc_orient", "2_mcn_pgs.rds")
    if (!dir.exists(output_dir)) stop("Directory ", output_dir, " does not exist.")
    if (file.exists(file.path(output_dir, "check_pc_number.log"))) {
        i <- 1
        while (file.exists(file.path(output_dir, paste0("check_pc_number.log.old", i)))) i <- i + 1
        invisible(file.rename(file.path(output_dir, "check_pc_number.log"),
                              file.path(output_dir, paste0("check_pc_number.log.old", i))))
    }
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
    if (!file.exists(info_path)) {
        logerror(paste0(info_path, " does not exist. Please run estimate_mcn_from_cel() with \"keep_tempfile = TRUE\" to generate the proper intermediate results for PC number check."))
        stop(info_path, "does not exist.")
    }
    if (max_pc_num <= 0 | (max_pc_num >= 1 & as.integer(max_pc_num) != max_pc_num)) {
        logerror("\"max_pc_num\" should be a fraction between 0 and 1 (proportion of post-QC sample size) or an positive integer (count of PCs to be used).")
        stop("Invalid \"max_pc_num\" argument.")
    }
    mito_lrr <- readRDS(mito_path)
    sample_qc_info <- readRDS(info_path)
    ind.pcincl <- sample_qc_info$used_in_PCA
    N_pcincl <- sum(ind.pcincl)
    if (is.null(correlated_pheno)) {
        if (N_pcincl >= 1000) {
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
        correlated_pheno <- correlated_pheno %>%
            select(Sample_Name, Phenotype) %>%
            drop_na()
        temp_common_samp <- sum(colnames(mito_lrr) %in% correlated_pheno$Sample_Name)
        if (temp_common_samp == 0) {
            logerror("No common sample between genotyping results and \"correlated_pheno\". Please check the sample names. (Sample names are defined as CEL file names. Did you forget the \".CEL\" surfix?)")
            stop("No common sample between genotyping results and \"correlated_pheno\".")
        }
        if (temp_common_samp < 4) {
            logerror("No enouch common sample between genotyping results and \"correlated_pheno\". Please check the sample names. (Sample names are defined as CEL file names. Did you forget the \".CEL\" surfix?)")
            stop("No enouch common sample between genotyping results and \"correlated_pheno\".")
        }
        if (temp_common_samp < 200) {
            logwarn(paste0("Only ", temp_common_samp, " common sample found between genotyping results and \"correlated_pheno\". The correlation may not have enough power."))
        }
    }
    loginfo("All preparation complete.")
    if (max_pc_num < 1) {
        numpc <- round(max_pc_num * N_pcincl)
        if (numpc == 0) numpc <- 1
    } else {
        numpc <- max_pc_num
    }
    if (max_pc_num < 1 & numpc < 15) {
        logwarn(paste0("Only ", numpc, "PC(s) were used (", N_pcincl, "*", max_pc_num, "). Consider set a higher \"max_pc_num\" to better evaluate PC numbers to be used."))
    }
    if (numpc > N_pcincl) {
        logwarn(paste0("Number of PCs used (", numpc, ") is larger than the post QC sample size (", N_pcincl, "). Will use ", N_pcincl, " PCs."))
        numpc <- N_pcincl
    }
    flag_calcPC <- TRUE
    if (file.exists(autopc_path)) {
        auto_pca_scores <- readRDS(autopc_path)
        if (ncol(auto_pca_scores) >= numpc) {
            loginfo("Autosome PCA file already exist. Skipping PCA...")
            auto_pca_scores <- auto_pca_scores[, 1:numpc]
            flag_calcPC <- FALSE
        } else {
            loginfo(paste0("Existing autosome PCA file doesn't have enough PCs (", ncol(auto_pca), "<", numpc, "). Recalculating ", numpc, " PC(s)..."))
        }
    } else {
        loginfo(paste0("Calculating ", numpc, " PC(s)..."))
    }
    if (flag_calcPC) {
        X_auto_lrr <- big_attach(auto_path)
        auto_pca <- big_randomSVD(X_auto_lrr, 
                                  fun.scaling = big_scale(), 
                                  ind.row = which(ind.pcincl), 
                                  k = numpc, 
                                  ncores = nb_cores())
        auto_pca_scores <- predict(auto_pca, X_auto_lrr)
        logdebug("Saving autosomal PC(s)...")
        saveRDS(auto_pca_scores, autopc_path)
    }
    # use only HQ samples
    auto_pca <- auto_pca_scores[ind.pcincl, ]
    mito_lrr <- mito_lrr[, ind.pcincl, with = F]
    loginfo("Calculating MCN estimates adjusted by different number of autosomal PCs...")
    registerDoParallel(detectCores())
    temp_res <- foreach(i = 1:numpc, .combine = c) %dopar% {
        mito_lrr_adjusted <- apply(mito_lrr, 1, function(X) {
            lm(X ~ auto_pca[, 1:i])$residuals
        })
        mito_pc1 <- as_FBM(mito_lrr_adjusted) %>%
            big_randomSVD(fun.scaling = big_scale(), 
                          k = 1, 
                          ncores = 1) %>%
            predict()
        mcn_raw <- data.frame(Sample_Name = rownames(mito_lrr_adjusted),
                              MCN = mito_pc1[, 1])
        temp_merge <- merge(mcn_raw, correlated_pheno)
        cor(temp_merge$Phenotype, temp_merge$MCN, method = "spearman")
    }
    stopImplicitCluster()
    plot(abs(temp_res), xlab = "Number of PC(s) included", ylab = paste0("Correlation between MCN estimates and ", pheno_name))
    loginfo("PC number check complete.")
}