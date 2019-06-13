## Assembly and SNP calling benchmarking based on HCMV mixed strain genome sequence samples

This repository contains the scripts and pipeline that reproduces the results of the HCMV benchmarking study. In this study we evaluated genome assemblers and variant callers on 6 in vitro generated, mixed strain HCMV sequence samples, each consisting of two lab strains in different abundance ratios. 

### Prerequirements

To reproduce the output, the following programs with specified versions need to be installed:

1. `Snakemake` (v5.3.0) for excuting the pipeline.
2. `R-3.5.1` and package `tidyverse`, `SomaticSignatures`, `MutationalPatterns` and `cowplot` for the visualization and mutation signature analysis.
3. `fastp` (v0.19.4), `BWA MEM` (v0.7.17), `SAMtools` (v1.9), `mummer` (v3.23), `picard` (v2.18.2), `csvtk` (v0.17.0).
4. `BCFtools` (v1.9), `VarScan` (v2.4.3), `Freebayes` (v1.2.0), `LoFreq` (v2.1.3.1).
5. `ABySS` (v2.1.4), `megahit` (v1.1.3) , `IDBA` (v1.1.3), `SPAdes` (v3.12.0), `Ray` (v2.3.1), `tadpole` (v37.99).
6. `metaquast` (v5.0.2) was used for evaluate the assembles.
7. `PEAR` (v0.9.6), `savage` (v0.4.0) for haplotype assembly. `Python` (v2.7) is required by `savage`.



First of all, `bioconda` must be installed on your system. You can follow these instruction [here](https://bioconda.github.io). Then, you can simply install all programs using:

```shell
conda env create -f config/conda_env.yaml
```

The above command line will create a conda environment named **hcmv_benchmark** and install all required programs listed above into this environment except for `savage`.

Since `savage` requires `Python` version 2.7, you need to create a distinct environment for it:

```shell
conda create -n savage savage=0.4.0 PEAR=0.9.6
```

Install `R` and `Python` packages:

```shell
# Activate environment
conda activate hcmv_benchmark
# Install Python packages via pip
pip install pandas
# Install R packages 
Rscript -e 'install.packages(c("tidyverse", "cowplot", "pheatmap", "VennDiagram", 
   "gridExtra", "reshape2", "NMF","ggdendro", "RColorBrewer"))'
```

To Install `R Bioconductor` packages `SomaticSignatures`, `MutationalPatterns` , we should type `R` to enter the `R` command window. In the `R`command window, you need to install core packages for `Bioconductor`:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

Then, the `SomaticSignatures`, `MutationalPatterns` can be installed by:

```R
BiocManager::install(c("SomaticSignatures", "MutationalPatterns", "SomaticCancerAlterations"))
```



After all this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:hzi-bifo/HCMV_benchmark.git
```

### Download the dataset and reference genomes
All sequencing data can be obtained from ENV with accession number: PRJEB32127. And the genome sequences of HCMV strains are included in this repo.
```shell
# Download the sequencing data
```


### Preapre the configuration file and sample list file

- Modify the `config/config.yaml`
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

```yaml
samplesDesc: config/sample.path.tsv
MerlinRef: ../ref/merlin/Merlin.BAC.fasta # path to Merlin genome
TB40ERef: ../ref/TB40/TB40E.GFP.fasta # path to TB40 genome
AD169Ref: ../ref/AD169/AD169.BAC.fasta # path to Ad169 genome
projectPath: ../hcmv_benchmark_output # path to the project aka the directory for the outputs

genomeDiffM: merlin_TB40E
genomeDiffA: AD169_TB40E

```

- Provide the sample list `config/sample_list.txt`. The list is a tab delimitaed text file, and echo row is one sample.
The format is as follows:
```tsv
sample	r1	r2
TA-0-1	../data/cleaned/cl_fq/TA-0-1.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-0-1.qc.nophix.r2.fq
TA-1-0	../data/cleaned/cl_fq/TA-1-0.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-1-0.qc.nophix.r2.fq
TA-1-10	../data/cleaned/cl_fq/TA-1-10.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-1-10.qc.nophix.r2.fq
TA-1-1	../data/cleaned/cl_fq/TA-1-1.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-1-1.qc.nophix.r2.fq
TA-1-50	../data/cleaned/cl_fq/TA-1-50.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-1-50.qc.nophix.r2.fq
TM-0-1	../data/cleaned/cl_fq/TM-0-1.qc.nophix.r1.fq	../data/cleaned/cl_fq/TM-0-1.qc.nophix.r2.fq
TM-1-0	../data/cleaned/cl_fq/TM-1-0.qc.nophix.r1.fq	../data/cleaned/cl_fq/TM-1-0.qc.nophix.r2.fq
TM-1-10	../data/cleaned/cl_fq/TM-1-10.qc.nophix.r1.fq	../data/cleaned/cl_fq/TM-1-10.qc.nophix.r2.fq
TM-1-1	../data/cleaned/cl_fq/TM-1-1.qc.nophix.r1.fq	../data/cleaned/cl_fq/TM-1-1.qc.nophix.r2.fq
TM-1-50	../data/cleaned/cl_fq/TM-1-50.qc.nophix.r1.fq	../data/cleaned/cl_fq/TM-1-50.qc.nophix.r2.fq
```
Please modify the paths to the sequencing files accordingly.

- Obtain the `CLC` SNP calling results from this repo 
`CLC` is not a freeware, so here we provide the output SNPs in this repo. First you need to create a project folder, for example with name `hcmv_benchmark_output` and then copy the files into the project directory:

```shell
mkdir -p ~/hcmv_benchmark_output/results/snp/callers/clc/
# put the CLC SNP calling results in the corresponding folder
cp -r HCMV_benchmark/data/clc ~/hcmv_benchmark_output/results/snp/callers/clc/
```

### Evaluate the assembly and haplotype reconstruction
```shell
conda activate hcmv_benchmark
snakemake -s evaluate_assembly.smk [-j <cores>]
```

### Test variant callers and analyze the mutation context of identified variants
To include the CLC result in the evaluation, please put the CLC variants calling results provided in this repo into your project directory:
```shell
mkdir -p <your project path>/results/SNP/callers
cp -r HCMV_benchmark/data/CLC <your project path>/results/SNP/callers
```
Evaluation
```shell
conda activate hcmv_benchmark
snakemake -s evaluate_snpcall.smk [-j <cores>]
```

### All in one:
```shell
python run_all_benchmarking.py
```
