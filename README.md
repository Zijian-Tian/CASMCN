# CASMCN

CASMCN is an R package for the estimation of relative mitochondrial copy number (MCN) using genotyping results (.CEL files) from [CAS Array](https://academic.oup.com/pcm/article/6/1/pbad002/7055961).

## System Requirements

In addition to the installation of the CASMCN package, [APT](https://www.thermofisher.com/cn/zh/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html) and [PennCNV](https://penncnv.openbioinformatics.org/en/latest/) are also needed for MCN estimation. The pipeline was tested on APT 2.11.4 and PennCNV 1.0.5, but should also work fine on other versions.

## Installation

CASMCN can be installed from github repository:

```r
# install.packages("remotes")
remotes::install_github("Zijian-Tian/CASMCN")
```

## Usage

A tutorial on estimating MCN using CASMCN is available [here](tutorial.md).

## Methods

This package was developed for estimating MCN for CAS Array, a customized Axiom array. The core processes of the package were developed with reference to two pipelines for MCN estimation, [MitoPipeline](http://genvisis.org/MitoPipeline/index.html) and [AutoMitoC](https://github.com/GMELab/AutoMitoC). 

In brief, this package estimate MCN in 3 steps. First, use APT to conduct quality controls on raw intensity data and extract the genotype calls and intensity signals. Then, use PennCNV to calculate log R ratio (LRR) indicating normalized intensity signals and conduct GC correction. Finally, use R to estimate MCN. High quality markers  selected by BLAST (exclude off-target probes) and PLINK (standard GWAS QC) are used to estimate MCN. Samples with excessive LRR fluctuation (LRR SD > 0.35) are excluded for low genotyping quality. PCA is conducted on LRR of selected autosomal markers and the top PCs capturing confounding factors are used to adjust the mitochondrial LRR. Another PCA on adjusted mitochondrial LRR generates the first PC as the measurement of relative MCN. The direction of the PC is determined by the correlation between MCN and an external phenotype correlated to MCN, such as age, sex, white blood cell counts or platelet count. If no external phenotype is provided, polygenic scores (PGS) of MCN calculated from the nuclear genotype is used to determine the direction of the first PC. The PGS model was generated by [LDpred2](https://github.com/privefl/bigsnpr) using [GWAS summary statistics from previous study](https://doi.org/10.1007/s00439-021-02394-w).

## Citation

[Zijian Tian, Fei Chen, Jing Wang, Benrui Wu, Jian Shao, Ziqing Liu, Li Zheng, You Wang, Tao Xu, Kaixin Zhou, CAS Array: design and assessment of a genotyping array for Chinese biobanking, Precision Clinical Medicine, Volume 6, Issue 1, March 2023, pbad002, doi:10.1093/pcmedi/pbad002](https://academic.oup.com/pcm/article/6/1/pbad002/7055961)

## References

MitoPipeline: Generating Mitochondrial copy number estimates from SNP array data in Genvisis. http://genvisis.org/MitoPipeline/index.html.

Chong M, Mohammadi-Shemirani P, Perrot N, Nelson W, Morton R, Narula S, et al. GWAS and ExWAS of blood mitochondrial DNA copy number identifies 71 loci and highlights a potential causal role in dementia. Elife. 2022;11.

Analysis Power Tools (APT). 2021. https://www.thermofisher.cn/cn/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html.

Wang K, Li M, Hadley D, Liu R, Glessner J, Grant SF, et al. PennCNV: an integrated hidden Markov model designed for high-resolution copy number variation detection in whole-genome SNP genotyping data. Genome Res. 2007;17(11):1665-74.

Diskin SJ, Li M, Hou C, Yang S, Glessner J, Hakonarson H, et al. Adjustment of genomic waves in signal intensities from whole-genome SNP genotyping platforms. Nucleic Acids Res. 2008;36(19):e126.

R Core Team. R: A Language and Environment for Statistical Computing. 2021. https://www.R-project.org/.

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: architecture and applications. BMC Bioinformatics. 2009;10(1):421.

Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7.

Privé F, Arbel J, Vilhjálmsson BJ. LDpred2: better, faster, stronger. Bioinformatics. 2020;36(22-23):5424-31.

Longchamps RJ, Yang SY, Castellani CA, Shi W, Lane J, Grove ML, et al. Genome-wide analysis of mitochondrial DNA copy number reveals loci implicated in nucleotide metabolism, platelet activation, and megakaryocyte proliferation. Human Genetics. 2022;141(1):127-46.
