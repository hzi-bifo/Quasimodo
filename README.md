## Assembly and SNP calling benchmarking based on HCMV mixed strain genome sequence samples

This repository contains the scripts and pipeline that reproduces the results of the HCMV benchmarking study. In this study we evaluated genome assemblers and variant callers on 6 in vitro generated, mixed strain HCMV sequence samples, each consisting of two lab strains in different abundance ratios. 

In this benchmarking: variants callers `BCFtools` (v1.9), `VarScan` (v2.4.3), `Freebayes` (v1.2.0), `LoFreq` (v2.1.3.1), `CLC Genomics Workbench` (v11.0.1) were evaluated. For the assembly benchmarking, `ABySS` (v2.1.4), `megahit` (v1.1.3) , `IDBA` (v1.1.3), `SPAdes` (v3.12.0), `Ray` (v2.3.1), `tadpole` (v37.99) were assessed. The haplotype reconstruction program `Savage` (v0.4.0) was also evaluated. 

### Prerequirements

To reproduce the output, you need to use `Bioconda`.

Please follow the instruction [here](https://bioconda.github.io) to install `Bioconda`. And then you need to install `snakemake`:

```shell
conda install snakemake=5.3.0
```

After this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:hzi-bifo/HCMV_benchmark.git
```

### Download the dataset and reference genomes
All sequencing data can be obtained from ENV with accession number: PRJEB32127. And the genome sequences of HCMV strains and Phix are included in this repo in the `ref` directory. Download the sequencing data into a directory, for example: `<your project path>/data/seqs/reads`:
```shell
# Create the dir for the reads
mkdir -p <your project path>/data/seqs/reads
# Download the sequencing data
cd HCMV_benchmark
wget -i data/PRJEB32127.txt -P <your project path>/data/seqs/reads
```

### Prepare the configuration file and sample list file

- Modify the `config/config.yaml`.
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

```yaml
samplesDesc: config/sample_list.tsv
MerlinRef: ref/Merlin.BAC.fa # path to Merlin genome
TB40ERef: ref/TB40E.GFP.fa # path to TB40 genome
AD169Ref: ref/AD169.BAC.fa # path to AD169 genome
PhixRef: ref/Phix.fa # path to Phix genome
projectPath: <your project path> # path to the project aka the directory for outputs
threads: 2 # number of cores to use
```

- Provide the sample list `config/sample_list.tsv`. The list is a tab delimited text file, and echo row is one sample.
The format is as follows:
```tsv
sample	r1	r2
TA-0-1	../HCMV_benchmark_output/data/seqs/reads/TA-0-1.qc.r1.fq.gz	../HCMV_benchmark_output/data/seqs/reads/TA-0-1.qc.r2.fq.gz
TA-1-0	../HCMV_benchmark_output/data/seqs/reads/TA-1-0.qc.r1.fq.gz	../HCMV_benchmark_output/data/seqs/reads/TA-1-0.qc.r2.fq.gz
TA-1-10	../HCMV_benchmark_output/data/seqs/reads/TA-1-10.qc.r1.fq.gz	../HCMV_benchmark_output/data/seqs/reads/TA-1-10.qc.r2.fq.gz
...
```
Please modify the paths to the sequencing files accordingly. In this example, the <your project path> is `../HCMV_benchmark_output` and the reads are in the `../HCMV_benchmark_output/data/seqs/reads`.


### Run the benchmarking

**All evaluation can be launched with `run_benchmark.py`**

#### Evaluate the assembly and haplotype reconstruction
```shell
python3 run_benchmark.py -t 10 assembly -c ~/miniconda3/envs
```
Please use -c parameter to specify the desired path to create the `conda` ENVs. The `conda` ENVs will be created under the path of the program by default. The program may take 1 hour to create the ENV used in this benchmarking for the first time.

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

```shell
python3 run_benchmark.py -t 10 snpcall -c ~/miniconda3/envs
```

#### All in one:
```shell
python run_benchmarking.py -t 10 all -c ~/miniconda3/envs
```
