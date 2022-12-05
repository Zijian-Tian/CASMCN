#' @title Estimate MCN from CEL files.
#' @description Generate MCN estimates from CEL files producted by CAS Array.
#' @param cel_files A character vector containing paths to all CEL files to be included in MCN estimation.
#' @param output_dir Path to the output directory.
#' @param apt_lib_dir Path to the library directory of CAS Array for APT. The directory should contain several XML for APT analysis and a CSV as array annotation.
#' @param gc5_file Path to the gzipped gc model file from PennCNV. The file could be found in gc_file/hg38.gc5Base.txt.gz in your local PennCNV installation or be downloaded at https://github.com/WGLab/PennCNV/blob/master/gc_file/hg38.gc5Base.txt.gz.
#' @param autosomal_markers (Optional) A character vector containing names of high-quality autosomal markers to be included in MCN estimation. Default is `hq_autosomal_markers` which contains a set of pre-defined high quality autosomal markers.
#' @param mitochondrial_markers (Optional) A character vector containing names of high-quality mitochondrial markers to be included in MCN estimation. Default is `hq_mitochondrial_markers` which contains a set of pre-defined high quality mitochondrial markers.
#' @param apt_exec_dir (Optional) Path to the directory containing APT executables. APT could be downloaded at https://www.thermofisher.com/cn/zh/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html. Default is `NULL` which means the APT executables will be searched from system `PATH`.
#' @param penncnv_exec_dir (Optional) Path to the directory containing PennCNV executables. PennCNV could be downloaded at https://github.com/WGLab/PennCNV/releases. Default is `NULL` which means the PennCNV executables will be searched from system `PATH`.
#' @param pc_used (Optional) A number indicating principal components to be used for confounders correction. Either a fraction between 0 and 1 indicating proportion of post-QC sample size or an positive integer indicating count of PCs to be used are accepted. Default is `0.01`.
#' @param keep_tempfile (Optional) Whether to keep the temporary files of intermediate results for debug purpose or further analysis? Default is `TRUE`.
#' @param ignore_existing_tempfiles (Optional) If detected existing temporary files in the output directory, whether to use the intermediate results or delete them before the analysis? Use intermediate results will save time significantly, but incomplete intermediate results could cause issues. Default is `FALSE`.
#' @param correlated_pheno (Optional) A data frame with columns "Sample_Name" and "Phenotype". The phenotype should be associated with MCN for the oriantation of the first PC in PCA of mitochondrial markers. Default is `NULL` which means the MCN PGS will be used for PC oriantation. Note that the association between MCN PGS and MCN estimates may not strong enough in small populations (e.g. less than a few hundred of samples) so an associated phenotype is highly recommended.
#' @param correlation_direction (Optional) A character indicating the direction of the association between the provided phenotype and MCN. Accepts "`+`" or "`-`" indicating positive or negative correlation between the phenotype and MCN. Default is "`+`".
#' @export
#' @import dplyr
#' @import logging
#' @import bigstatsr
#' @importFrom data.table fread setDTthreads
#' @importFrom HardyWeinberg HWExactStats
#' @importFrom stats cor.test lm sd predict
#' @importFrom utils read.csv read.table write.table
#' @importFrom tidyr drop_na
estimate_mcn_from_cel <- function(cel_files, output_dir, apt_lib_dir, gc5_file, autosomal_markers = hq_autosomal_markers, mitochondrial_markers = hq_mitochondrial_markers, apt_exec_dir = NULL, penncnv_exec_dir = NULL, pc_used = 0.01, keep_tempfile = TRUE, ignore_existing_tempfiles = FALSE, correlated_pheno = NULL, correlation_direction = '+') {
    if (!dir.exists(output_dir)) dir.create(output_dir)
    if (file.exists(file.path(output_dir, "estimate_MCN_from_CEL.log"))) {
        i <- 1
        while (file.exists(file.path(output_dir, paste0("estimate_MCN_from_CEL.log.old", i)))) i <- i + 1
        invisible(file.rename(file.path(output_dir, "estimate_MCN_from_CEL.log"),
                              file.path(output_dir, paste0("estimate_MCN_from_CEL.log.old", i))))
    }
    setDTthreads(0)
    logReset()
    addHandler(writeToConsole)
    addHandler(writeToFile, file = file.path(output_dir, "estimate_MCN_from_CEL.log"))
    logger <- getLogger()
    logger$setLevel('DEBUG')
    setLevel('INFO', logger$handlers$writeToConsole)
    setLevel('DEBUG', logger$handlers$writeToFile)
    loginfo("Starts MCN estimation from CEL files.")
    loginfo(paste0("Log file created at ", file.path(output_dir, "estimate_MCN_from_CEL.log"), "."))
    loginfo("Checking arguments...")
    if (pc_used <= 0 | (pc_used >= 1 & as.integer(pc_used) != pc_used)) {
        logerror("\"pc_used\" should be a fraction between 0 and 1 (proportion of post-QC sample size) or an positive integer (count of PCs to be used).")
        stop("Invalid \"pc_used\" argument.")
    }
    if (!(correlation_direction %in% c('+', '-'))) {
        logerror("\"correlation_direction\" should be '+' or '-'.")
        stop("Invalid \"correlation_direction\" argument.")
    }
    if (is.null(correlated_pheno)) {
        if (length(cel_files) >= 1000) {
            loginfo("No correlated phenotype provided. Will use MCN PRS for PC oriantation.")
        } else {
            logwarn("No correlated phenotype provided. Will use MCN PRS for PC oriantation. Note that this may not powerful enough in small sample size (e.g. <1000).")
        }
    } else {
        if (!all(c("Sample_Name", "Phenotype") %in% colnames(correlated_pheno))) {
            logerror("\"correlated_pheno\" should be either NULL or a data frame with \"Sample_Name\" and \"Phenotype\" columns.")
            stop("Invalid \"correlated_pheno\" argument.")
        }
        correlated_pheno <- correlated_pheno %>%
            select(Sample_Name, Phenotype) %>%
            drop_na()
        temp_common_samp <- sum(basename(cel_files) %in% correlated_pheno$Sample_Name)
        if (temp_common_samp == 0) {
            logerror("No common sample between \"cel_files\" and \"correlated_pheno\". Please check the sample names. (Sample names are defined as CEL file names. Did you forget the \".CEL\" surfix?)")
            stop("No common sample between \"cel_files\" and \"correlated_pheno\".")
        }
        if (temp_common_samp < 4) {
            logerror("No enouch common sample between \"cel_files\" and \"correlated_pheno\". Please check the sample names. (Sample names are defined as CEL file names. Did you forget the \".CEL\" surfix?)")
            stop("No enouch common sample between \"cel_files\" and \"correlated_pheno\".")
        }
        if (temp_common_samp < 200) {
            logwarn(paste0("Only ", temp_common_samp, " common sample found between \"cel_files\" and \"correlated_pheno\". The correlation may not have enough power."))
        }
    }
    logdebug("Checking if every file in CEL list exists...")
    for (cel_file in cel_files) {
        if (!file.exists(cel_file)) {
            logerror(paste(cel_file, "in CEL list is not found."))
            stop(paste(cel_file, "in CEL list is not found."))
        }
    }
    if (length(cel_files) < 1000) {
        logwarn("MCN estimates may not powerful enough in small populations.")
    }
    logdebug("Checking if APT library directory exists...")
    if (!dir.exists(apt_lib_dir)) {
        logerror(paste(apt_lib_dir, "does not exist."))
        stop(paste(apt_lib_dir, "does not exist."))
    }
    logdebug("Checking if every APT library file needed exists...")
    temp_libs <- list.files(apt_lib_dir)
    temp_dqc <- grepl("\\.apt-geno-qc\\.AxiomQC1\\.xml$", temp_libs)
    if (sum(temp_dqc) == 0) {
        logerror(paste0("DQC library file not found at ", apt_lib_dir, "."))
        stop(paste0("DQC library file not found at ", apt_lib_dir, "."))
    } else if (sum(temp_dqc) > 1) {
        logwarn(paste0("Multiple DQC library file found at ", apt_lib_dir, ". The first one will be used."))
        lib_dqc <- temp_libs[temp_dqc]
        lib_dqc <- lib_dqc[1]
    } else {
        lib_dqc <- temp_libs[temp_dqc]
    }
    lib_dqc <- file.path(apt_lib_dir, lib_dqc)
    logdebug(paste0("Using DQC library file at ", lib_dqc, "."))
    temp_crqc <- grepl("\\_96orMore\\_Step1", temp_libs) & 
                 grepl("\\.apt-genotype-axiom\\.AxiomGT1\\.apt2\\.xml$", temp_libs)
    if (sum(temp_crqc) == 0) {
        logerror(paste0("Call Rate QC library file not found at ", apt_lib_dir, "."))
        stop(paste0("Call Rate QC library file not found at ", apt_lib_dir, "."))
    } else if (sum(temp_crqc) > 1) {
        logwarn(paste0("Multiple Call Rate QC library file found at ", apt_lib_dir, ". The first one will be used."))
        lib_crqc <- temp_libs[temp_crqc]
        lib_crqc <- lib_crqc[1]
    } else {
        lib_crqc <- temp_libs[temp_crqc]
    }
    lib_crqc <- file.path(apt_lib_dir, lib_crqc)
    logdebug(paste0("Using Call Rate QC library file at ", lib_crqc, "."))
    temp_gt <- grepl("\\_96orMore\\_Step2", temp_libs) &
               grepl("\\.apt-genotype-axiom\\.", temp_libs) &
               grepl("\\.AxiomGT1\\.apt2\\.xml$", temp_libs)
    if (sum(temp_gt) == 0) {
        logerror(paste0("Genotyping library file not found at ", apt_lib_dir, "."))
        stop(paste0("Genotyping library file not found at ", apt_lib_dir, "."))
    } else if (sum(temp_gt) > 1) {
        logwarn(paste0("Multiple genotyping library file found at ", apt_lib_dir, ". The first one will be used."))
        lib_gt <- temp_libs[temp_gt]
        lib_gt <- lib_gt[1]
    } else {
        lib_gt <- temp_libs[temp_gt]
    }
    lib_gt <- file.path(apt_lib_dir, lib_gt)
    logdebug(paste0("Using genotyping library file at ", lib_gt, "."))
    temp_anno <- grepl("\\.annot\\.csv$", temp_libs)
    anno_zipped <- FALSE
    if (sum(temp_anno) == 0) {
        temp_anno <- grepl("\\.annot\\.csv\\.zip$", temp_libs)
        if (sum(temp_anno) == 0) {
            logerror(paste0("Annotation library file not found at ", apt_lib_dir, "."))
            stop(paste0("Annotation library file not found at ", apt_lib_dir, "."))
        } else if (sum(temp_anno) > 1) {
            anno_zipped <- TRUE
            logwarn(paste0("Multiple annotation library file found at ", apt_lib_dir, ". The first one will be used."))
            lib_anno <- temp_libs[temp_anno]
            lib_anno <- lib_anno[1]
        } else {
            anno_zipped <- TRUE
            lib_anno <- temp_libs[temp_anno]
        }
    } else if (sum(temp_anno) > 1) {
        logwarn(paste0("Multiple annotation library file found at ", apt_lib_dir, ". The first one will be used."))
        lib_anno <- temp_libs[temp_anno]
        lib_anno <- lib_anno[1]
    } else {
        lib_anno <- temp_libs[temp_anno]
    }
    lib_anno <- file.path(apt_lib_dir, lib_anno)
    logdebug(paste0("Using annotation library file at ", lib_anno, "."))
    logdebug("All APT library files exist.")
    logdebug("Checking gc5Base file...")
    if (!file.exists(gc5_file)) {
        logerror(paste0("Gc5Base file not found at ", apt_lib_dir, "."))
        stop(paste0("Gc5Base file not found at ", apt_lib_dir, "."))
    }
    logdebug("Checking APT excutables...")
    if (!is.null(apt_exec_dir)) {
        if (!dir.exists(apt_exec_dir)) {
            logerror(paste0("APT excutable directory ", apt_exec_dir, " does not exist."))
            stop(paste0("APT excutable directory ", apt_exec_dir, " does not exist."))
        } else {
            if (!file.exists(file.path(apt_exec_dir, "apt-geno-qc"))) {
                logerror(paste0("apt-geno-qc does not exist at APT excutable directory ", apt_exec_dir, "."))
                stop(paste0("apt-geno-qc does not exist at APT excutable directory ", apt_exec_dir, "."))
            }
            if (!file.exists(file.path(apt_exec_dir, "apt-genotype-axiom"))) {
                logerror(paste0("apt-genotype-axiom does not exist at APT excutable directory ", apt_exec_dir, "."))
                stop(paste0("apt-genotype-axiom does not exist at APT excutable directory ", apt_exec_dir, "."))
            }
            if (!file.exists(file.path(apt_exec_dir, "ps-extract"))) {
                logerror(paste0("ps-extract does not exist at APT excutable directory ", apt_exec_dir, "."))
                stop(paste0("ps-extract does not exist at APT excutable directory ", apt_exec_dir, "."))
            }
        }
    } else {
        for (i in c("apt-geno-qc", "apt-genotype-axiom", "ps-extract")) {
            if (system(paste0("command -v ", i, " > /dev/null 2>&1")) != 0) {
                logerror(paste0("Can't find ", i, " in system PATH. Please check APT is installed and added to PATH or set the APT excutable directory explicitly."))
                stop(paste0("Can't find ", i, " in system PATH."))
            }
        }
    }
    logdebug("Checking PennCNV excutables...")
    if (!is.null(penncnv_exec_dir)) {
        if (!dir.exists(penncnv_exec_dir)) {
            logerror(paste0("PennCNV excutable directory ", penncnv_exec_dir, " does not exist."))
            stop(paste0("PennCNV excutable directory ", penncnv_exec_dir, " does not exist."))
        } else {
            if (!file.exists(file.path(penncnv_exec_dir, "cal_gc_snp.pl"))) {
                logerror(paste0("cal_gc_snp.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "."))
                stop(paste0("cal_gc_snp.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "."))
            }
            if (!file.exists(file.path(penncnv_exec_dir, "kcolumn.pl"))) {
                logerror(paste0("kcolumn.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "."))
                stop(paste0("kcolumn.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "."))
            }
            if (!file.exists(file.path(penncnv_exec_dir, "genomic_wave.pl"))) {
                logerror(paste0("genomic_wave.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "."))
                stop(paste0("genomic_wave.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "."))
            }
            if (!file.exists(file.path(penncnv_exec_dir, "affy/bin/generate_affy_geno_cluster.pl"))) {
                logerror(paste0("generate_affy_geno_cluster.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "/affy/bin."))
                stop(paste0("generate_affy_geno_cluster.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "/affy/bin."))
            }
            if (!file.exists(file.path(penncnv_exec_dir, "affy/bin/normalize_affy_geno_cluster.pl"))) {
                logerror(paste0("normalize_affy_geno_cluster.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "/affy/bin."))
                stop(paste0("normalize_affy_geno_cluster.pl does not exist at PennCNV excutable directory ", penncnv_exec_dir, "/affy/bin."))
            }
        }
    } else {
        for (i in c("cal_gc_snp.pl", "kcolumn.pl", "genomic_wave.pl", "generate_affy_geno_cluster.pl", "normalize_affy_geno_cluster.pl")) {
            if (system(paste0("command -v ", i, " > /dev/null 2>&1")) != 0) {
                logerror(paste0("Can't find ", i, " in system PATH. Please check PennCNV is installed and added to PATH or set the PennCNV excutable directory explicitly. Note that scripts processing affy arrays are nested in affy/bin/ of PennCNV excutable directory so additional entry of PATH is needed."))
                stop(paste0("Can't find ", i, " in system PATH."))
            }
        }
    }
    loginfo("Arguments checking complete.")

    temp_dir <- file.path(output_dir, "tempdir")
    if (dir.exists(temp_dir)) {
        if (ignore_existing_tempfiles) {
            loginfo("Deleting existing tempfiles...")
            unlink(temp_dir, recursive = T)
            dir.create(temp_dir)
        } else {
            logwarn("Intermediate results from existing tempfiles will be used. If the arguments have changed or the intermediate results are incomplete, this may cause issues. You can delete incomplete intermediate results manually to re-run substeps or set \"ignore_existing_tempfiles\" to TRUE to re-run the whole pipeline.")
        }
    } else {
        logdebug("No existing tempfiles found.")
        dir.create(temp_dir)
    }

    loginfo("All preparation complete.")

    loginfo("STEP1: Call genotype using APT")
    geno_dir <- file.path(temp_dir, "1_APT_genotyping")
    dir.create(geno_dir, showWarnings = F)
    geno_tempdir <- file.path(geno_dir, "tempdir")
    dir.create(geno_tempdir, showWarnings = F)
    logdebug("Writing CEL list...")
    cel_full <- file.path(geno_tempdir, "cel_full.txt")
    write.table(data.frame(cel_files = cel_files), cel_full, quote = F, row.names = F, col.names = T)
    loginfo("STEP1.1: DQC with APT")
    dqc_dir <- file.path(geno_dir, "1_dqc")
    dir.create(dqc_dir, showWarnings = F)
    dqc_out <- file.path(dqc_dir, "CASMCN_apt-geno-qc.report.txt")
    dqc_log <- file.path(dqc_dir, "CASMCN_apt-geno-qc_dqc.log")
    if (file.exists(dqc_out)) {
        loginfo("DQC output file already exists. Skipping STEP1.1...")
    } else {
        temp_exec <- "apt-geno-qc"
        if (!is.null(apt_exec_dir)) temp_exec <- file.path(apt_exec_dir, temp_exec)
        temp_cmd <- paste(temp_exec,
            "--analysis-files-path", apt_lib_dir,
            "--xml-file", lib_dqc,
            "--out-dir", dqc_dir,
            "--cel-files", cel_full,
            "--out-file", dqc_out,
            "--log-file", dqc_log)
        logdebug(paste0("Call command: ", temp_cmd))
        system(temp_cmd, ignore.stdout = T, ignore.stderr = T)
        if (!file.exists(dqc_out)) {
            logerror(paste0("DQC output file not found at ", dqc_out, ". Please check APT log at ", dqc_log, "."))
            stop(paste0("DQC output file not found at ", dqc_out, "."))
        }
    }
    dqc_cel <- read.table(dqc_out, header = T, sep = "\t", comment.char = "#") %>% 
        select(cel_files, axiom_dishqc_DQC) %>%
        filter(axiom_dishqc_DQC >= 0.82) %>%
        select(cel_files)
    logdebug("DQC output file found. Writing CEL list passed DQC...")
    cel_dqc <- file.path(geno_tempdir, "cel_pass_dqc.txt")
    data.frame(cel_files = cel_files) %>%
        filter(basename(cel_files) %in% dqc_cel$cel_files) %>%
        write.table(cel_dqc, quote = F, row.names = F, col.names = T)
    loginfo(paste0(nrow(dqc_cel), " CEL files passed DQC."))
    loginfo("STEP1.2: Calling Rate QC with APT")
    crqc_dir <- file.path(geno_dir, "2_crqc")
    dir.create(crqc_dir, showWarnings = F)
    crqc_out <- file.path(crqc_dir, "AxiomGT1.report.txt")
    crqc_log <- file.path(crqc_dir, "CASMCN_apt-genotype-axiom_crqc.log")
    if (file.exists(crqc_out)) {
        loginfo("Call Rate QC output file already exists. Skipping STEP1.2...")
    } else {
        temp_exec <- "apt-genotype-axiom"
        if (!is.null(apt_exec_dir)) temp_exec <- file.path(apt_exec_dir, temp_exec)
        temp_cmd <- paste(temp_exec,
            "--table-output", "true",
            "--analysis-files-path", apt_lib_dir,
            "--arg-file", lib_crqc,
            "--cel-files", cel_dqc,
            "--out-dir", crqc_dir,
            "--log-file", crqc_log)
        logdebug(paste0("Calling command: ", temp_cmd))
        system(temp_cmd, ignore.stdout = T, ignore.stderr = T)
        if (!file.exists(crqc_out)) {
            logerror(paste0("Call Rate QC output file not found at ", crqc_out, ". Please check APT log at ", crqc_log, "."))
            stop(paste0("Call Rate QC output file not found at ", crqc_out, "."))
        }
    }
    crqc_cel <- read.table(crqc_out, header = T, sep = "\t", comment.char = "#") %>% 
        select(cel_files, call_rate) %>%
        filter(call_rate >= 97.0) %>%
        select(cel_files)
    logdebug("Call Rate QC output file found. Writing CEL list passed Call Rate QC...")
    cel_crqc <- file.path(geno_tempdir, "cel_pass_crqc.txt")
    data.frame(cel_files = cel_files) %>%
        filter(basename(cel_files) %in% crqc_cel$cel_files) %>%
        write.table(cel_crqc, quote = F, row.names = F, col.names = T)
    loginfo(paste0(nrow(crqc_cel), " CEL files passed Call Rate QC."))
    loginfo("STEP1.3: Genotyping with APT")
    gt_dir <- file.path(geno_dir, "3_genotyping")
    dir.create(gt_dir, showWarnings = F)
    if (all(file.exists(file.path(gt_dir, "AxiomGT1.calls.txt")),
            file.exists(file.path(gt_dir, "AxiomGT1.confidences.txt")),
            file.exists(file.path(gt_dir, "AxiomGT1.summary.txt")))) {
        loginfo("Genotyping output files already exist. Skipping STEP1.3...")
    } else {
        temp_exec <- "apt-genotype-axiom"
        if (!is.null(apt_exec_dir)) temp_exec <- file.path(apt_exec_dir, temp_exec)
        gt_log <- file.path(gt_dir, "CASMCN_apt-genotype-axiom_genotyping.log")
        temp_cmd <- paste(temp_exec,
            "--analysis-files-path", apt_lib_dir,
            "--arg-file", lib_gt,
            "--cel-files", cel_crqc,
            "--out-dir", gt_dir,
            "--log-file", gt_log,
            "--allele-summaries", "true")
        logdebug(paste0("Calling command: ", temp_cmd))
        system(temp_cmd, ignore.stdout = T, ignore.stderr = T)
        for (i in c("AxiomGT1.calls.txt", "AxiomGT1.confidences.txt", "AxiomGT1.summary.txt")) {
            if (!file.exists(file.path(gt_dir, i))) {
                logerror(paste0("Genotyping output files not found at ", file.path(gt_dir, i), ". Please check APT log at ", gt_log, "."))
                stop(paste0("Genotyping output files not found at ", file.path(gt_dir, i), "."))
            }
        }
    }
    loginfo(paste0("Genotyping complete for ", nrow(crqc_cel), " samples."))
    loginfo("STEP1.4: Extract genotyping results with APT")
    ext_dir <- file.path(geno_dir, "4_extract")
    dir.create(ext_dir, showWarnings = F)
    if (all(file.exists(file.path(ext_dir, "extract_calls.txt")),
            file.exists(file.path(ext_dir, "extract_confidences.txt")),
            file.exists(file.path(ext_dir, "extract_summary.txt")))) {
        loginfo("Extracted genotype files already exist. Skipping STEP1.4...")
    } else {
        pid_full <- file.path(geno_tempdir, "all_probesets.txt")
        data.frame(probeset_id = c(autosomal_markers, mitochondrial_markers)) %>%
            write.table(pid_full, quote = F, row.names = F, col.names = T)
        ext_log <- file.path(ext_dir, "CASArray_extract.log")
        temp_exec <- "ps-extract"
        if (!is.null(apt_exec_dir)) temp_exec <- file.path(apt_exec_dir, temp_exec)
        temp_cmd <- paste(temp_exec,
            "--call-file", file.path(gt_dir, "AxiomGT1.calls.txt"),
            "--confidence-file", file.path(gt_dir, "AxiomGT1.confidences.txt"),
            "--summary-file", file.path(gt_dir, "AxiomGT1.summary.txt"),
            "--pid-file", pid_full,
            "--output-dir", ext_dir,
            "--log-file", ext_log)
        logdebug(paste0("Calling command: ", temp_cmd))
        system(temp_cmd, ignore.stdout = T, ignore.stderr = T)
        for (i in c("extract_calls.txt", "extract_confidences.txt", "extract_summary.txt")) {
            if (!file.exists(file.path(ext_dir, i))) {
                logerror(paste0("Extracted genotyping output file not found at ", file.path(ext_dir, i), "."))
                stop(paste0("Extracted genotyping output file not found at ", file.path(ext_dir, i), "."))
            }
        }
    }
    loginfo(paste0("Genotype extracted for ", nrow(crqc_cel), " samples."))

    loginfo("STEP2: Call LRR using PennCNV")
    lrr_dir <- file.path(temp_dir, "2_PennCNV_LRR")
    dir.create(lrr_dir, showWarnings = F)
    lrr_tempdir <- file.path(lrr_dir, "tempdir")
    dir.create(lrr_tempdir, showWarnings = F)
    gcm_gc5 <- file.path(lrr_tempdir, "gc5Base.txt.sorted")
    penn_mapfile <- file.path(lrr_tempdir, "CASArray_mapfile.dat")
    penn_mapfileax <- file.path(lrr_tempdir, "CASArray_mapfileAX.dat")
    penn_snpfile <- file.path(lrr_tempdir, "CASArray_snpfile.txt")
    if (all(file.exists(gcm_gc5), file.exists(penn_mapfile), file.exists(penn_mapfileax), file.exists(penn_snpfile))) {
        loginfo("Processed library files already exist. Skipping preprocessing...")
    } else {
        loginfo("Prepare files...")
        logdebug("Processing GC annotation file...")
        conn <- gzfile(gc5_file, 'rt')
        conn %>%
            read.table(header = F) %>%
            arrange(V2, V3) %>%
            write.table(gcm_gc5, sep = "\t", quote = F, col.names = F, row.names = F)
        close(conn)
        logdebug("Processing array annotation file...")
        if (anno_zipped) {
            conn <- unz(lib_anno, gsub("\\.zip", "", basename(lib_anno)))
            anno_full <- read.csv(conn, comment.char = "#")
            close(conn)
        } else {
            anno_full <- read.csv(lib_anno, comment.char = "#")
        }
        anno_proc <- anno_full %>%
            filter(Probe.Set.ID %in% c(autosomal_markers, mitochondrial_markers)) %>%
            select(Probe.Set.ID, Chromosome, Physical.Position) %>%
            rename(Name = Probe.Set.ID, Chr = Chromosome, Position = Physical.Position) %>%
            mutate(Chr = gsub("MT", "M", Chr))
        write.table(select(anno_proc, Name, Chr), penn_mapfile, sep = '\t', quote = F, row.names = F, col.names = F)
        write.table(anno_proc, penn_mapfileax, sep = '\t', quote = F, row.names = F, col.names = F)
        write.table(anno_proc, penn_snpfile, sep = '\t', quote = F, row.names = F, col.names = T)
    }
    loginfo("STEP2.1: Generate GC-model with PennCNV")
    gcm_dir <- file.path(lrr_dir, "1_gc_model")
    dir.create(gcm_dir, showWarnings = F)
    gcm_out <- file.path(gcm_dir, "CASArray.gcmodel")
    gcm_log <- file.path(gcm_dir, "CASArray.gcmodel.log")
    if (file.exists(gcm_out)) {
        loginfo("GC-model output file already exists. Skipping STEP2.1...")
    } else {
        temp_exec <- "cal_gc_snp.pl"
        if (!is.null(penncnv_exec_dir)) temp_exec <- file.path(penncnv_exec_dir, temp_exec)
        temp_cmd <- paste(temp_exec, gcm_gc5, penn_snpfile, "-output", gcm_out, ">", gcm_log, "2>&1")
        logdebug(paste0("Call command: ", temp_cmd))
        system(temp_cmd)
        if (!file.exists(gcm_out)) {
            logerror(paste0("PennCNV GC calculation output file not found at ", gcm_out, ". Please check log at ", gcm_log, "."))
            stop(paste0("PennCNV GC calculation output file not found at ", gcm_out, "."))
        }
    }
    loginfo("STEP2.2: Genotype clustering with PennCNV")
    cluster_dir <- file.path(lrr_dir, "2_clustering")
    dir.create(cluster_dir, showWarnings = F)
    cluster_out <- file.path(cluster_dir, "CASArray.genocluster")
    cluster_log <- file.path(cluster_dir, "CASArray.genocluster.log")
    if (file.exists(cluster_out)) {
        loginfo("Genotype clustering file already exists. Skipping STEP2.2...")
    } else {
        temp_exec <- "generate_affy_geno_cluster.pl"
        if (!is.null(penncnv_exec_dir)) temp_exec <- file.path(penncnv_exec_dir, "affy", "bin", temp_exec)
        temp_cmd <- paste(temp_exec, 
                        file.path(ext_dir, "extract_calls.txt"),
                        file.path(ext_dir, "extract_confidences.txt"),
                        file.path(ext_dir, "extract_summary.txt"),
                        "-nopower2",
                        "-locfile", penn_mapfile,
                        "-out", cluster_out,
                        ">", cluster_log, "2>&1")
        logdebug(paste0("Call command: ", temp_cmd))
        system(temp_cmd)
        if (!file.exists(cluster_out)) {
            logerror(paste0("PennCNV genotype clustering output file not found at ", cluster_out, ". Please check log at ", cluster_log, "."))
            stop(paste0("PennCNV genotype clustering output file not found at ", cluster_out, "."))
        }
    }
    loginfo("STEP2.3: LRR calculation with PennCNV")
    lrrcalc_dir <- file.path(lrr_dir, "3_LRR_calc")
    dir.create(lrrcalc_dir, showWarnings = F)
    lrrcalc_out <- file.path(lrrcalc_dir, "CASArray_lrr_baf.txt")
    lrrcalc_log <- file.path(lrrcalc_dir, "CASArray_lrr_baf.log")
    if (file.exists(lrrcalc_out)) {
        loginfo("LRR calculation output file already exists. Skipping STEP2.3...")
    } else {
        temp_exec <- "normalize_affy_geno_cluster.pl"
        if (!is.null(penncnv_exec_dir)) temp_exec <- file.path(penncnv_exec_dir, "affy", "bin", temp_exec)
        temp_cmd <- paste(temp_exec, 
                        cluster_out,
                        file.path(ext_dir, "extract_summary.txt"),
                        "-nopower2",
                        "-locfile", penn_mapfileax,
                        "-out", lrrcalc_out,
                        ">", lrrcalc_log, "2>&1")
        logdebug(paste0("Call command: ", temp_cmd))
        system(temp_cmd)
        if (!file.exists(lrrcalc_out)) {
            logerror(paste0("PennCNV LRR calculation output file not found at ", lrrcalc_out, ". Please check log at ", lrrcalc_log, "."))
            stop(paste0("PennCNV LRR calculation output file not found at ", lrrcalc_out, "."))
        }
    }
    loginfo("STEP2.4: Split LRR&BAF with PennCNV")
    lrrsplit_dir <- file.path(lrr_dir, "4_LRR_split")
    dir.create(lrrsplit_dir, showWarnings = F)
    lrrsplit_split_dir <- file.path(lrrsplit_dir, "splitted")
    dir.create(lrrsplit_split_dir, showWarnings = F)
    lrrsplit_out <- file.path(lrrsplit_dir, "CASArray_lrr_split.filenames")
    lrrsplit_log <- file.path(lrrsplit_dir, "CASArray_lrr_split.log")
    lrrsplit_out_LRR <- file.path(lrrsplit_dir, "CASArray_lrr_split_LRR.filenames")
    if (file.exists(lrrsplit_out_LRR)) {
        loginfo("Split LRR output file already exists. Skipping STEP2.4...")
    } else {
        temp_exec <- "kcolumn.pl"
        if (!is.null(penncnv_exec_dir)) temp_exec <- file.path(penncnv_exec_dir, temp_exec)
        temp_cmd <- paste(temp_exec, 
                        lrrcalc_out,
                        "split", "1",
                        "--tab",
                        "--head", "3",
                        "--name_by_header",
                        "--output", file.path(lrrsplit_split_dir, "LRR_split"),
                        "--filenameout", lrrsplit_out,
                        ">", lrrsplit_log, "2>&1")
        logdebug(paste0("Call command: ", temp_cmd))
        system(temp_cmd)
        if (!file.exists(lrrsplit_out)) {
            logerror(paste0("PennCNV LRR kcolumn split output file not found at ", lrrsplit_out, ". Please check log at ", lrrsplit_log, "."))
            stop(paste0("PennCNV LRR kcolumn split output file not found at ", lrrsplit_out, "."))
        }
        read.table(lrrsplit_out, header = F) %>% 
            filter(!grepl("\\_R1$", V1)) %>% 
            write.table(lrrsplit_out_LRR, quote = F, col.names = F, row.names = F)
    }
    loginfo("STEP2.5: GC correction with PennCNV")
    gccorr_dir <- file.path(lrr_dir, "5_GC_correction")
    dir.create(gccorr_dir, showWarnings = F)
    gccorr_log <- file.path(gccorr_dir, "CASArray_GC_correction.log")
    gccorr_out <- file.path(gccorr_dir, "CASArray_GC_adjusted.filelist")
    if (file.exists(gccorr_out)) {
        loginfo("GC correction output file already exists. Skipping STEP2.5...")
        gccorr_files <- read.table(gccorr_out, header = F)$V1
    } else {
        temp_exec <- "genomic_wave.pl"
        if (!is.null(penncnv_exec_dir)) temp_exec <- file.path(penncnv_exec_dir, temp_exec)
        temp_cmd <- paste(temp_exec, 
                        "-adjust",
                        "-gcmodel", gcm_out,
                        "--listfile", lrrsplit_out_LRR,
                        ">", gccorr_log, "2>&1")
        logdebug(paste0("Call command: ", temp_cmd))
        system(temp_cmd)
        gccorr_files <- list.files(lrrsplit_split_dir, full.names = T, pattern = "\\.adjusted$")
        if (length(gccorr_files) == 0) {
            logerror(paste0("PennCNV GC correction output file not found. Please check log at ", gccorr_log, "."))
            stop(paste0("PennCNV GC correction output file not found."))
        }
        data.frame(files = gccorr_files) %>%
            write.table(gccorr_out, quote = F, row.names = F, col.names = F)
    }

    loginfo("STEP3: Estimate MCN using R")
    mcn_dir <- file.path(temp_dir, "3_R_MCN")
    dir.create(mcn_dir, showWarnings = F)
    if (all(file.exists(file.path(mcn_dir, "3_sample_QC.rds")),
            file.exists(file.path(mcn_dir, "3_autosomal_LRR.rds")),
            file.exists(file.path(mcn_dir, "3_autosomal_LRR.bk")),
            file.exists(file.path(mcn_dir, "3_mitochondrial_LRR.rds")))) {
        loginfo("Extracted LRR files already exist. Skipping STEP3.1-3.3...")
        sample_qc_info <- readRDS(file.path(mcn_dir, "3_sample_QC.rds"))
        X_auto_lrr <- big_attach(file.path(mcn_dir, "3_autosomal_LRR.rds"))
        mito_lrr <- readRDS(file.path(mcn_dir, "3_mitochondrial_LRR.rds"))
    } else {
        loginfo("STEP3.1: Read LRR files")
        all_lrr_rowinfo <- fread(gccorr_files[1], header = TRUE, sep = "\t", select = 1:3)
        loginfo(paste0("Reading ", sum(all_lrr_rowinfo$Chr != 'M'), " autosomal markers and ", sum(all_lrr_rowinfo$Chr == 'M'), " mitochondrial markers from LRR files..."))
        all_lrr <- do.call(cbind, lapply(gccorr_files, function(x) fread(x,  header = T, sep = '\t', select = 4)))
        all_lrr_colnames <- colnames(all_lrr)
        all_lrr_colnames <- sub(".Log R Ratio", "", all_lrr_colnames)
        colnames(all_lrr) <- all_lrr_colnames
        loginfo(paste0("Read LRR of ", length(all_lrr_colnames), " samples from LRR files."))
        loginfo("STEP3.2: QC of selected markers")
        n_crqc <- length(gccorr_files)
        snp_info <- fread(cluster_out, header = T, select = c("probeset_id", "count_aa", "count_ab", "count_bb"))
        snp_info_mt <- filter(snp_info, probeset_id %in% mitochondrial_markers)
        ind.cr.mt <- snp_info_mt$count_aa + snp_info_mt$count_ab + snp_info_mt$count_bb < 0.95 * n_crqc
        logdebug(paste0(sum(ind.cr.mt), " mitochondrial markers excluded for low call rate (<0.95)."))
        mt_marker_used <- snp_info_mt$probeset_id[!ind.cr.mt]
        loginfo(paste0(length(mt_marker_used), " mitochondrial markers were included for MCN estimation."))
        snp_info_auto <- filter(snp_info, probeset_id %in% autosomal_markers)
        ind.cr <- snp_info_auto$count_aa + snp_info_auto$count_ab + snp_info_auto$count_bb < 0.95 * n_crqc
        snp_info_auto <- snp_info_auto[!ind.cr, ]
        logdebug(paste0(sum(ind.cr), " autosomal markers excluded for low call rate (<0.95)."))
        baf <- (2 * snp_info_auto$count_ab + snp_info_auto$count_bb) / (snp_info_auto$count_aa + snp_info_auto$count_ab + snp_info_auto$count_bb) / 2
        ind.maf <- baf < 0.01 | baf > 0.99
        snp_info_auto <- snp_info_auto[!ind.maf, ]
        logdebug(paste0(sum(ind.maf), " autosomal markers excluded for low allele frequency (<0.01)."))
        ind.hwe <- HWExactStats(snp_info_auto[, 2:4]) < 1e-6
        snp_info_auto <- snp_info_auto[!ind.hwe, ]
        logdebug(paste0(sum(ind.hwe), " autosomal markers excluded for low HWE p-value (<1e-6)."))
        auto_marker_used <- snp_info_auto$probeset_id
        loginfo(paste0(length(auto_marker_used), " autosomal markers were included for MCN estimation."))
        write.table(data.frame(x = auto_marker_used), file.path(mcn_dir, "2_post_QC_markers_autosome.txt"), quote = F, row.names = F, col.names = F)
        write.table(data.frame(x = mt_marker_used), file.path(mcn_dir, "2_post_QC_markers_mitochondrial.txt"), quote = F, row.names = F, col.names = F)
        loginfo("STEP3.3: QC of sample LRR")
        ind.auto <- which(all_lrr_rowinfo$Name %in% auto_marker_used)
        auto_lrr <- all_lrr[ind.auto, ]
        lrr_sd <- apply(auto_lrr, 2, sd)
        ind.lrrsd <- lrr_sd > 0.35
        sample_qc_info <- data.frame(sample_name = all_lrr_colnames, lrr_sd = lrr_sd, used_in_PCA = !ind.lrrsd)
        write.table(sample_qc_info, file.path(mcn_dir, "3_sample_QC.txt"), quote = F, sep = '\t', col.names = T, row.names = F)
        saveRDS(sample_qc_info, file.path(mcn_dir, "3_sample_QC.rds"))
        logdebug(paste0(sum(ind.lrrsd), " samples excluded for high LRR SD (>0.35)."))
        loginfo(paste0(sum(!ind.lrrsd), " samples were included for MCN estimation."))
        ind.mito <- which(all_lrr_rowinfo$Name %in% mt_marker_used)
        mito_lrr <- all_lrr[ind.mito, ]
        logdebug("Saving LRR...")
        if (file.exists(file.path(mcn_dir, "3_autosomal_LRR.bk"))) {
            logdebug("Deleting old FBM...")
            unlink(file.path(mcn_dir, "3_autosomal_LRR.bk"))
        }
        X_auto_lrr <- as_FBM(t(auto_lrr), backingfile = file.path(mcn_dir, "3_autosomal_LRR"))
        saveRDS(X_auto_lrr, file.path(mcn_dir, "3_autosomal_LRR.rds"))
        saveRDS(mito_lrr, file.path(mcn_dir, "3_mitochondrial_LRR.rds"))
        rm(all_lrr)
    }
    loginfo("STEP3.4: PCA of autosomal LRR")
    ind.pcincl <- sample_qc_info$used_in_PCA
    N_pcincl <- sum(ind.pcincl)
    if (pc_used < 1) {
        numpc <- round(pc_used * N_pcincl)
        if (numpc == 0) numpc <- 1
    } else {
        numpc <- pc_used
    }
    if (pc_used < 1 & numpc < 15) {
        logwarn(paste0("Only ", numpc, "PC(s) were used (", N_pcincl, "*", pc_used, "). Consider set \"pc_used\" manually to better capture confounders."))
    }
    if (numpc > N_pcincl) {
        logwarn(paste0("Number of PCs used (", numpc, ") is larger than the post QC sample size (", N_pcincl, "). Will use ", N_pcincl, " PCs."))
        numpc <- N_pcincl
    }
    flag_calcPC <- TRUE
    if (file.exists(file.path(mcn_dir, "4_autosomal_PCA.rds"))) {
        auto_pca_scores <- readRDS(file.path(mcn_dir, "4_autosomal_PCA.rds"))
        if (ncol(auto_pca_scores) >= numpc) {
            loginfo("Autosome PCA file already exist. Skipping STEP3.4...")
            auto_pca_scores <- auto_pca_scores[, 1:numpc]
            flag_calcPC <- FALSE
        } else {
            loginfo(paste0("Existing autosome PCA file doesn't have enough PCs (", ncol(auto_pca_scores), "<", numpc, "). Recalculating ", numpc, " PC(s)..."))
        }
    } else {
        loginfo(paste0("Calculating ", numpc, " PC(s)..."))
    }
    if (flag_calcPC) {
        auto_pca <- big_randomSVD(X_auto_lrr, 
                                  fun.scaling = big_scale(), 
                                  ind.row = which(ind.pcincl), 
                                  k = numpc, 
                                  ncores = nb_cores())
        auto_pca_scores <- predict(auto_pca, X_auto_lrr)
        logdebug("Saving autosomal PC(s)...")
        saveRDS(auto_pca_scores, file.path(mcn_dir, "4_autosomal_PCA.rds"))
    }
    loginfo("STEP3.5: Adjust mitochondrial LRR by autosomal PCs")
    mito_lrr_adjusted <- apply(mito_lrr, 1, function(X) {
        lm(X ~ auto_pca_scores)$residuals
    })
    loginfo("STEP3.6: PCA of mitochondrial LRR")
    X_mito_lrr_adjusted <- as_FBM(mito_lrr_adjusted)
    mito_pca <- big_randomSVD(X_mito_lrr_adjusted, 
                              fun.scaling = big_scale(), 
                              ind.row = which(ind.pcincl), 
                              k = 1, 
                              ncores = nb_cores())
    mito_pc1 <- predict(mito_pca, X_mito_lrr_adjusted)
    mcn_raw <- data.frame(Sample_Name = rownames(mito_lrr_adjusted),
                          MCN = mito_pc1[, 1],
                          used_in_PCA = ind.pcincl)
    loginfo("STEP3.7: Orientation of mitochondrial PC")
    if (!is.null(correlated_pheno)) {
        if (correlation_direction == '+') {
            loginfo("Orientating PC with phenotype positively correlated with MCN...")
        } else {
            loginfo("Orientating PC with phenotype negatively correlated with MCN...")
        }
        temp_merge <- merge(mcn_raw[ind.pcincl, ], correlated_pheno, by = "Sample_Name")
        temp_cor <- cor.test(temp_merge$MCN, temp_merge$Phenotype, method = "spearman")
        if (temp_cor$p.value > 0.05) logwarn(paste0("Correlation between MCN and reference phenotype provided is not statistically significant (p = ", signif(temp_cor$p.value, 2), " > 0.05). Consider using other phenotypes or increasing the sample size."))
        if ((correlation_direction == '+' & temp_cor$estimate > 0) |
            (correlation_direction == '-' & temp_cor$estimate < 0)) {
            mcn_final <- mcn_raw
        } else {
            mcn_final <- mcn_raw %>% mutate(MCN = -1 * MCN)
        }
        final_output <- file.path(output_dir, "CASMCN_result.txt")
        write.table(mcn_final, final_output, quote = F, sep = '\t', row.names = F, col.names = T)
        loginfo(paste0("Final MCN estimate saved at ", final_output, "."))
        temp_merge <- merge(mcn_final, correlated_pheno, by = "Sample_Name")
        temp_cor <- cor.test(temp_merge$MCN, temp_merge$Phenotype, method = "spearman")
        loginfo(paste0("Correlation between MCN and reference phenotype: Spearman's Rho = ", signif(temp_cor$estimate, 2), ", p = ", signif(temp_cor$p.value, 2), "."))
    } else {
        loginfo("Will orientate PC with MCN PGS. Preparing for PGS calculation...")
        prs_dir <- file.path(mcn_dir, "7_pc_orient")
        dir.create(prs_dir, showWarnings = F)
        loginfo("Extracting genotyping results...")
        prs_ext_dir <- file.path(prs_dir, "1_extract")
        dir.create(prs_ext_dir, showWarnings = F)
        prs_mcn <- file.path(prs_dir, "2_mcn_pgs.rds")
        if (file.exists(prs_mcn)) {
            loginfo("MCN polygenic score file already exists. Skipping PGS calculation...")
            prs_df <- readRDS(prs_mcn)
        } else {
            prs_ext_geno <- file.path(prs_ext_dir, "extract_calls.txt")
            if (file.exists(prs_ext_geno)) {
                loginfo("Extracted genotype files already exist. Skipping extraction...")
            } else {
                prs_pid <- file.path(prs_dir, "prs_probesets.txt")
                mcnpgs_betas %>%
                    select(probeid) %>%
                    rename(probeset_id = probeid) %>%
                    write.table(prs_pid, quote = F, row.names = F, col.names = T)
                prs_ext_log <- file.path(prs_ext_dir, "CASArray_extract.log")
                temp_exec <- "ps-extract"
                if (!is.null(apt_exec_dir)) temp_exec <- file.path(apt_exec_dir, temp_exec)
                temp_cmd <- paste(temp_exec,
                    "--call-file", file.path(gt_dir, "AxiomGT1.calls.txt"),
                    "--pid-file", prs_pid,
                    "--output-dir", prs_ext_dir,
                    "--log-file", prs_ext_log)
                logdebug(paste0("Calling command: ", temp_cmd))
                system(temp_cmd, ignore.stdout = T, ignore.stderr = T)
                if (!file.exists(prs_ext_geno)) {
                    logerror(paste0("Extracted genotyping output file not found at ", prs_ext_geno, "."))
                    stop(paste0("Extracted genotyping output file not found at ", prs_ext_geno, "."))
                }
            }
            temp_geno <- fread(prs_ext_geno, sep = "\t", header = T, skip = "probeset_id")
            loginfo("Extracted genotype obtained. Performing QC for genotypes...")
            temp_geno %>%
                select(-probeset_id) %>%
                as.matrix -> geno_mat
            rownames(geno_mat) <- temp_geno$probeset_id
            rm(temp_geno)
            geno_mat[geno_mat < 0] <- NA
            ind.samp <- apply(geno_mat, 2, function(X) sum(is.na(X))) / nrow(geno_mat) > 0.05
            logdebug(paste0(sum(ind.samp), " samples excluded from PGS calculation due to low call rate (<0.05)."))
            geno_mat <- geno_mat[, !ind.samp]
            ind.cr <- apply(geno_mat, 1, function(X) sum(is.na(X))) / ncol(geno_mat) > 0.05
            logdebug(paste0(sum(ind.cr), " variants excluded from PGS calculation due to low call rate (<0.05)."))
            geno_mat <- geno_mat[!ind.cr, ]
            mafs <- rowMeans(geno_mat, na.rm = T) / 2
            ind.maf <- (mafs < 0.01) | (mafs > 0.99)
            logdebug(paste0(sum(ind.maf), " variants excluded from PGS calculation due to low MAF (<0.01)."))
            geno_mat <- geno_mat[!ind.maf, ]
            var_count <- cbind(apply(geno_mat, 1, function(X) sum(X == 0, na.rm = T)),
                            apply(geno_mat, 1, function(X) sum(X == 1, na.rm = T)),
                            apply(geno_mat, 1, function(X) sum(X == 2, na.rm = T)))
            ind.hwe <- HWExactStats(var_count) < 1e-6
            logdebug(paste0(sum(ind.hwe), " variants excluded from PGS calculation due to low HWE p-value (<1e-6)."))
            geno_mat <- geno_mat[!ind.hwe, ]
            ind.na <- which(is.na(geno_mat), arr.ind = T)
            geno_mat[ind.na] <- rowMeans(geno_mat, na.rm = T)[ind.na[, "row"]]
            loginfo(paste0("Calculating PGS of MCN with ", ncol(geno_mat), " samples and ", nrow(geno_mat), " variants..."))
            beta_vec <- merge(data.frame(probeid = rownames(geno_mat)), mcnpgs_betas)$beta
            prs_res <- beta_vec %*% geno_mat
            prs_df <- data.frame(Sample_Name = colnames(prs_res), PRS = prs_res[, ])
            saveRDS(prs_df, prs_mcn)
        }
        loginfo("Orientating PC with MCN PGS...")
        temp_merge <- merge(prs_df, mcn_raw[ind.pcincl, ])
        temp_cor <- cor.test(temp_merge$MCN, temp_merge$PRS, method = "spearman")
        if (temp_cor$p.value > 0.05) logwarn(paste0("Correlation between MCN and PGS is not statistically significant (p = ", signif(temp_cor$p.value, 2), " > 0.05). Consider using other phenotypes or increasing the sample size."))
        if (temp_cor$estimate > 0) {
            mcn_final <- mcn_raw
        } else {
            mcn_final <- mcn_raw %>% mutate(MCN = -1 * MCN)
        }
        final_output <- file.path(output_dir, "CASMCN_result.txt")
        write.table(mcn_final, final_output, quote = F, sep = '\t', row.names = F, col.names = T)
        loginfo(paste0("Final MCN estimate saved at ", final_output, "."))
        temp_merge <- merge(prs_df, mcn_final)
        temp_cor <- cor.test(temp_merge$MCN, temp_merge$PRS, method = "spearman")
        loginfo(paste0("Correlation between MCN and reference PGS: Spearman's Rho = ", signif(temp_cor$estimate, 2), ", p = ", signif(temp_cor$p.value, 2), "."))
    }
    loginfo("MCN estimation from CEL files complete.")
    if (!keep_tempfile) {
        loginfo("Deleting temporary files...")
        unlink(temp_dir, recursive = T)
    }
}
