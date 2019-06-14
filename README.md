## Assembly and SNP calling benchmarking based on HCMV mixed strain genome sequence samples

This repository contains the scripts and pipeline that reproduces the results of the HCMV benchmarking study. In this study we evaluated genome assemblers and variant callers on 6 in vitro generated, mixed strain HCMV sequence samples, each consisting of two lab strains in different abundance ratios. 

In this benchmarking: variants callers `BCFtools` (v1.9), `VarScan` (v2.4.3), `Freebayes` (v1.2.0), `LoFreq` (v2.1.3.1), `CLC Genomics Workbench 11.0.1` were evaluated. For the assembly benchmarking, `ABySS` (v2.1.4), `megahit` (v1.1.3) , `IDBA` (v1.1.3), `SPAdes` (v3.12.0), `Ray` (v2.3.1), `tadpole` (v37.99) were assessed. The haplotype reconstruction program `Savage` (v0.4.0) was also evaluated. 

### Prerequirements

To reproduce the output, you need to use `Bioconda`.

Please follow the instruction [here](https://bioconda.github.io) to install `Bioconda`. After this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:hzi-bifo/HCMV_benchmark.git
```

### Download the dataset and reference genomes
All sequencing data can be obtained from ENV with accession number: PRJEB32127. And the genome sequences of HCMV strains are included in this repo in the `ref` directory.
```shell
# Download the sequencing data
wget -i data/PRJEB32127.txt -P <your project path>/data/seqs/reads
```

### Preapre the configuration file and sample list file

- Modify the `config/config.yaml`.
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

```yaml
samplesDesc: config/sample.path.tsv
MerlinRef: ../ref/merlin/Merlin.BAC.fasta # path to Merlin genome
TB40ERef: ../ref/TB40/TB40E.GFP.fasta # path to TB40 genome
AD169Ref: ../ref/AD169/AD169.BAC.fasta # path to AD169 genome
projectPath: <your project path> # path to the project aka the directory for outputs
threads: 2
```

- Provide the sample list `config/sample_list.tsv`. The list is a tab delimitaed text file, and echo row is one sample.
The format is as follows:
```tsv
sample	r1	r2
TA-0-1	../data/cleaned/cl_fq/TA-0-1.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-0-1.qc.nophix.r2.fq
TA-1-0	../data/cleaned/cl_fq/TA-1-0.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-1-0.qc.nophix.r2.fq
TA-1-10	../data/cleaned/cl_fq/TA-1-10.qc.nophix.r1.fq	../data/cleaned/cl_fq/TA-1-10.qc.nophix.r2.fq
...
```
Please modify the paths to the sequencing files accordingly. The reads are in the `<your project path>/data/seqs/reads`

- Obtain the `CLC` SNP calling results from this repo. 
`CLC` is not a freeware, so here we provide the output SNPs in this repo. First you need to create a project folder (which is defined in the config file as `projectPath`), for example with name `hcmv_benchmark_output` in your home directory and then copy the files into the project directory:

```shell
# enter the program directory
cd HCMV_benchmark
# create the directory
mkdir -p ~/hcmv_benchmark_output/results/snp/callers
# put the CLC SNP calling results in the corresponding folder
cp -r data/clc ~/hcmv_benchmark_output/results/snp/callers/
```

### Run the benchmarking

**All evaluation can be launched with `run_benchmark.py`**

#### Evaluate the assembly and haplotype reconstruction
```shell
python3 run_benchmark.py -t 10 assembly -c ~miniconda3/envs
```
The parameters for `run_benchmark.py`:
```shell
usage: run_benchmark.py [-h] [-d] [-t THREADS] [-c CONDA_PREFIX]
                        {all,snpcall,assembly}

positional arguments:
  {all,snpcall,assembly}
                        the evaluation to run

optional arguments:
  -h, --help            show this help message and exit
  -d, --dryrun          Print the output details without run the pipeline
  -t THREADS, --threads THREADS
                        The number of threads to use, default: 2
  -c CONDA_PREFIX, --conda_prefix CONDA_PREFIX
                        The prefix of conda ENV which tells the program where to create the conda ENV 
                        (default: in the working directory)
```

#### Test variant callers and analyze the mutation context of identified variants
To include the CLC result in the evaluation, please put the CLC variants calling results provided in this repo into your project directory:
```shell
mkdir -p <your project path>/results/SNP/callers
cp -r HCMV_benchmark/data/CLC <your project path>/results/SNP/callers
```
Evaluation
```shell
python3 run_benchmark.py -t 10 snpcall -c ~miniconda3/envs
```

#### All in one:
```shell
python run_benchmarking.py -t 10 all
```
