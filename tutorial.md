# MCN Estimation Tutorial

## Overview

This tutorial describes how to esitmate relative mitochondrial copy number (MCN) using genotyping results (.CEL files) from CAS Array by CASMCN package.

## Input Files

### CEL Files

Genotyping with CAS Array produces a CEL file for each sample. We use a character vector to locate all CEL files included in MCN estimation:

```r
cel_list <- c("/path/to/sample1.CEL",
              "/path/to/sample2.CEL",
              "/path/to/sample3.CEL",
              ...)
```

### Library Files

APT needs a set of library file to process CEL files. Library files of CAS Array could be obtained from manufacturer. The extracted library files should contain multiple .xml parameter files, a .annot.csv annotation file and several other files essential for APT processing. 

Most library files needed by PennCNV will be generaged from the annotation file of the CAS Array, but the gc model file need to be specified additionally. Since the genome build of the library files of CAS Array is hg38, the gc model file needs to be the hg38 version. You could find the file in gc_file/hg38.gc5Base.txt.gz in your local PennCNV installation, or download the file from [the github repository of PennCNV](https://github.com/WGLab/PennCNV/blob/master/gc_file/hg38.gc5Base.txt.gz).

```r
apt_lib_dir <- "/path/to/APT_library_directory"
gc5_file <- "/path/to/hg38.gc5Base.txt.gz"
```

## Basic Usage

The minimal requirements to run MCN estimation is the input files above. To start the analysis, use the function `estimate_mcn_from_cel()`:

```r
library(CASMCN)
estimate_mcn_from_cel(cel_list,
                      output_dir = "/path/to/output_directory",
                      apt_lib_dir,
                      gc5_file)
```

This function will generate relative MCN estimates for the samples and save the result to "CASMCN_result.txt" in your `output_dir`. A log file "estimate_MCN_from_CEL.log" will also be generated for tracking the progress since the analysis may take up to a few hundred of hours when applied to larger sample size. 

In addition to the basic options, **we recommend to provide a phenotype associated to MCN and choose an optimal number of PCs used (see below)** for the PC oriantation step, especially when the sample size is small (e.g. <1000).

## Output Format

The final result of MCN estimation will be written to "CASMCN_result.txt" in the output directory. It is a tab-delimited text file with two columns "Sample_Name" and "MCN". The sample name is the file name of each CEL file. MCN is the estimated relative mitochondrial copy number. **Note that the estimated MCN is not the absolute copy number of mitochondria. The relative MCN estimates could not be compared between different populations.**

## Optional Arguments

### Associated Phenotype

An external phenotype associated with MCN is needed to oriantate the PC derived from mitochondrial intensity signal. The correlated phenotype should be provided as a data frame with two columns "Sample_Name" and "Phenotype":

```r
head(correlated_pheno)
#   Sample_Name   Phenotype
# 1     sample1  0.87591389
# 2     sample2 -1.60092431
# 3     sample3  0.01470713
# 4     sample4  1.00549282
# 5     sample5 -0.09102592
# 6     sample6  0.35284596
```

The phenotype should have known correlation with MCN. Commonly used traits includes white blood cell counts (negative correlation), red blood cel counts (negative correlation), platelet count (positive correlation), age (negative correlation) and sex (female > male). If you have available MCN estimates for part of the samples, that also could be used in this step (as a phenotype positively associted with MCN). The direction of the association need to be clarified by `correlation_direction`:

```r
correlation_direction <- "+"  # phenotype positively correlated with MCN
# OR
correlation_direction <- "-"  # phenotype negatively correlated with MCN
```

If no external phenotype is provided (`correlated_pheno = NULL`), PC oriantation will be done with a polygenic score (PGS) of MCN derived from a [GWAS of MCN](https://doi.org/10.1007/s00439-021-02394-w) using [LDpred2](https://doi.org/10.1093/bioinformatics/btaa1029). Since the SNP heritability of MCN is not high enough, PGS may not powerful enough for PC oriantation in small population (e.g. <1000). **We recommend to provide external phenotype whenever it is possible.**

### Autosomal Principal Components Included

We use principal components (PCs) of intensity signal from autosomal markers to capture confounding factors including batch effects and DNA concentrations. If the number of PCs is not enough to capture the confounders, the result is unreliable. Empirically, 1% of the sample size is enough for MCN estimation in large populations. For MCN estimation in smaller populations, number of PCs may need to be set manually.

To control the number of PCs included, you can set `pc_used` in two ways: By proportion of the sample size passed the previous quality control (a fraction between 0 and 1); Or by the exact number of PCs included (a positive integer).

```r
pc_used <- 0.02  # By proportion of sample size
# OR
pc_used <- 20    # By exact number of PCs
```

To determine how many PCs should be used, **we recommend to check the MCN adjusted with different number of PCs by `plot_pc_adjustment()` function**. This function takes the `output_dir` of `estimate_mcn_from_cel()` as the input and use its intermediate results to plot correlation between estimated MCN and external phenotype (or MCN PGS if no phenotype is provided) against number of PCs used for adjustment. A typical plot is as follows:

![pc_adj_plot][pc_adj_plot]

The association between estimated MCN and external phenotype become stable after adjusted for ~100 PCs, giving a suggestion for the PC number selection. 

The `pc_used` parameter of the `estimate_mcn_from_cel()` function can be adjusted safely while using existing temporary files from previous runs: With full intermediate results, you can just adjust the `pc_used` and re-run `estimate_mcn_from_cel()` to generate a new MCN estimate. This should be fast if the new `pc_used` is smaller than the previous one since all necessery information could be loaded from the intermediate results. If the specified number of PCs used is larger than the cached PCs, another PCA had to be done. Since PCA is a time-consuming process, set a larger number of PCs used in the first run is recommended. That will cache enough information for further analysis, including determine the optimal number of PCs by `plot_pc_adjustment()` and estimate MCN using `estimate_mcn_from_cel()` with different settings of `pc_used`.

### High Quality Markers

A set of high quality markers for MCN estimate are specified by `autosomal_markers` and `mitochondrial_markers`. By default, a set of pre-defined high quality markers is used: 

```r
str(hq_autosomal_markers)
# chr [1:48107] "AX-32416607" "AX-40357457" "AX-32609917" "AX-32680447" ...
str(hq_mitochondrial_markers)
# chr [1:292] "AX-106719108" "AX-11414493" "AX-11414494" "AX-11414495" ...
```

To generate the pre-defined high quality markers, we applied the following filters: Exclude markers with off-target probes (more than one alignment with percentage of identical matches > 80%) for autosomal and mitochondrial markers; For autosomal markers, additional quality controls including call rate > 95%, HWE p-value > 10<sup>-6</sup>, LD r<sup>2</sup> < 0.3 and maximum spacing were done. 

### Temporary Files

This pipeline will generate a series of temporary files since multiple external softwares are included. All temporary files will be put in the "tempdir\" directory in the output directory. If `keep_tempfile` is set to `FALSE`, then all temporary files will be deleted after the final MCN estimates are generated. 

The pipeline could utilize legacy temporary files generate by previous runs. Unless `ignore_existing_tempfiles` is set to `TRUE`, the pipeline will skip the substeps as long as it detected the output of that substep. If that intermediate result is incomplete, this may cause issues. You can delete the output of the last substep and re-run the pipeline to solve the issue. If `ignore_existing_tempfiles` is set to `TRUE`, all temporary files will be deleted before the pipeline starts.

### External Software Directories

By default, the executables of APT and PennCNV will be searched from system `PATH`. You can also specify the main directory of APT and PennCNV installation explicitly by `apt_exec_dir` and `penncnv_exec_dir` arguments.

[pc_adj_plot]:data:image/png;base64,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