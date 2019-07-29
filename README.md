## Assembly and SNP calling benchmarking based on HCMV mixed strain genome sequence samples

This repository contains the scripts and pipeline that reproduces the results of the HCMV benchmarking study. In this study we evaluated genome assemblers and variant callers on 6 in vitro generated, mixed strain HCMV sequence samples, each consisting of two lab strains in different abundance ratios. 

In this benchmarking: variants callers `BCFtools` (v1.9), `VarScan` (v2.4.3), `Freebayes` (v1.2.0), `LoFreq` (v2.1.3.1), `CLC Genomics Workbench` (v11.0.1) were evaluated. For the assembly benchmarking, `ABySS` (v2.1.4), `megahit` (v1.1.3) , `IDBA` (v1.1.3), `SPAdes` (v3.12.0), `Ray` (v2.3.1), `tadpole` (v37.99) were assessed. The haplotype reconstruction program `Savage` (v0.4.0) was also evaluated. 

### Prerequirements

To reproduce the output, you need to use `Bioconda`.

Please follow the instruction [here](https://bioconda.github.io) to install `Bioconda`. And then you need to install `snakemake` and Python package `click`:

```shell
conda install snakemake=5.3.0
conda install click=7.0
```

After this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:hzi-bifo/HCMV_benchmark.git
```

### Download the dataset and prepare the sample list file (optional, only in the case you specify `--slow` option)
If you want to run the whole analyses on the reads you should download the reads. All sequencing data can be obtained from ENV with accession number: PRJEB32127. Download the sequencing data into a directory, for example: `<your project path>/data/seqs/reads`:
```shell
# Create the dir for the reads
mkdir -p <your project path>/data/seqs/reads
# Download the sequencing data
cd HCMV_benchmark
wget -i data/PRJEB32127.txt -P <your project path>/data/seqs/reads
```

- Provide the sample list `config/sample_list.tsv`. The list is a tab delimited text file, and echo row is one sample. The format is as follows:
```tsv
sample	r1	r2
TA-0-1	../HCMV_benchmark_output/data/seqs/reads/TA-0-1.qc.r1.fq.gz	../HCMV_benchmark_output/data/seqs/reads/TA-0-1.qc.r2.fq.gz
TA-1-0	../HCMV_benchmark_output/data/seqs/reads/TA-1-0.qc.r1.fq.gz	../HCMV_benchmark_output/data/seqs/reads/TA-1-0.qc.r2.fq.gz
TA-1-10	../HCMV_benchmark_output/data/seqs/reads/TA-1-10.qc.r1.fq.gz	../HCMV_benchmark_output/data/seqs/reads/TA-1-10.qc.r2.fq.gz
...
```
Please modify the paths to the sequencing files which you have downloaded accordingly. In this example, the `<your project path>` is `../HCMV_benchmark_output` and the reads are in the `../HCMV_benchmark_output/data/seqs/reads`.


#### ! Due to the high computational and time cost, by default this program do not run the whole benchmark for HCMV dataset from scratch (based on reads), instead it benchmarks the SNP call and assembly based on the VCF files and scaffolds provided within this program under `data` directory. 

### Adapt the configuration file
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

The genome sequences of HCMV strains and Phix are included in this repo in the `ref` directory. 

- Modify the `config/config.yaml`:

```yaml
samplesDesc: config/sample_list.tsv # leave it empty if you do not use --slow option
MerlinRef: ref/Merlin.BAC.fa # path to Merlin genome
TB40ERef: ref/TB40E.GFP.fa # path to TB40 genome
AD169Ref: ref/AD169.BAC.fa # path to AD169 genome
PhixRef: ref/Phix.fa # path to Phix genome
outpath: <your output path> # the directory for outputs
threads: 2 # number of cores to use
runOnReads: false # Run the whole analyses on reads. Controlled by the `--slow` option
```


### Run the benchmarking

**All evaluation can be launched with `run_benchmark.py`**

The parameters for `run_benchmark.py`:
```shell
Usage: run_benchmark.py [OPTIONS] COMMAND [ARGS]...

Options:
  --version      Print the version.
  --help         Show this message and exit.

Commands:
  hcmv     Benchmarking for HCMV dataset
  snpeval  SNP calling benchmark for customized dataset
  asmeval  Assembly benchmark for customized dataset
```

This program consists of three subcommands: `hcmv`, `snpeval`, `asmeval`. The first one is used for the benchmarking on our HCMV datasets. And the other two are for the SNP call and assembly evaluation on customized datasets.

The argumentrs and options in the `hcmv` command:
```shell
Usage: run_benchmark.py hcmv [OPTIONS]

  Benchmarking for HCMV dataset

Options:
  -o, --outpath PATH              The directory where to put the results and
                                  figures. The path can be specified either in
                                  the CLI as argument or in the config file.
  -c, --conda_prefix PATH         The prefix of conda ENV. [default: in the
                                  working directory].
  -t, --threads INTEGER           The number of threads to use.  [default: 2]
  -d, --dryrun                    Print the details without run the pipeline.
                                  [default: False]
  -e, --evaluation [all|snpcall|assembly]
                                  The evaluation to run.  [required]
  -s, --slow                      Run the evaluation based on reads, which is
                                  very slow. By default, the evaluation will
                                  be based on the VCF and contig files
                                  provided within this software. If this
                                  parameter is on, this software will run all
                                  the analyses to generate outputs based on
                                  reads for benchmarking which is very time
                                  consuming.  [default: False]
  --help                          Show this message and exit.
```


#### Evaluate the assembly and haplotype reconstruction
```shell
python3 run_benchmark.py hcmv -e assembly -t 10 -c ~/miniconda3/envs
```
Please use -c parameter to specify the desired path to create the `conda` ENVs. The `conda` ENVs will be created under the path of the program by default. The program may take 1 hour to create the ENV used in this benchmarking for the first time.

If you expect to the benchmarking based on the reads, you need to specify the `--slow` or `-s` option which allows you to generate the assembly results from reads.


#### Assess variant callers and analyze the mutation context of identified variants

```shell
python3 run_benchmark.py hcmv -e snpcall -t 10 -c ~/miniconda3/envs
```
If you wish to the benchmarking based on the reads, you need to specify the `--slow` or `-s` option which allows you to generate the variant calling results from reads.


#### All in one:
```shell
python run_benchmarking.py hcmv -e all -t 10 -c ~/miniconda3/envs
```

### The output structure

```
└── results
    ├── asssembly
    │   ├── abyss
    │   ├── idba
    │   ├── megahit
    │   ├── metaspades
    │   ├── ray
    │   ├── savage
    │   ├── spades
    │   └── tadpole
    ├── final_figures
    ├── final_tables
    ├── metaquast
    │   ├── summary_for_figure
    │   ├── TA
    │   └── TM
    └── snp
        ├── callers
        └── nucmer
```

#### For full run based on reads (-s option on)
```
├── data
│   └── seqs
│       ├── bam
│       ├── cl_bam
│       ├── cl_fq
│       ├── pear_merge
│       ├── pileup
│       └── reads
├── reports
│   ├── benchmarks
│   ├── bwa
│   └── picard
└── results
    ├── asssembly
    │   ├── abyss
    │   ├── idba
    │   ├── megahit
    │   ├── metaspades
    │   ├── ray
    │   ├── savage
    │   ├── spades
    │   └── tadpole
    ├── final_figures
    │   └── refined
    ├── final_tables
    ├── metaquast
    │   ├── summary_for_figure
    │   ├── TA
    │   └── TM
    └── snp
        ├── callers
        └── nucmer
```

### Run the benchmark on your own datasets

#### Evaluate the assembly or haplotype reconstruction

The arguments and options of `asmeval` command:
```shell
Usage: run_benchmark.py asmeval [OPTIONS]

  Assembly benchmark for customized dataset

Options:
  -o, --outpath PATH       The directory where to put the results and figures.
                           The path can be specified either in the CLI as
                           argument or in the config file.
  -c, --conda_prefix PATH  The prefix of conda ENV. [default: in the working
                           directory].
  -t, --threads INTEGER    The number of threads to use.  [default: 2]
  -d, --dryrun             Print the details without run the pipeline.
                           [default: False]
  -s, --scaffolds TEXT     Comma-separated list of scaffold files. Please
                           quote the whole parameter if there is any white
                           space the file names. The files can be specified
                           either in the CLI as argument or in the config
                           file.
  -r, --refs TEXT          Comma-separated list of reference genome files.
                           Please quote the whole parameter if there is any
                           white space the file names. The files can be
                           specified either in the CLI as argument or in the
                           config file.
  --help                   Show this message and exit.
```

- Run the benchmarking
```shell
python3 run_benchmark.py asmeval -t 10 -c ~/miniconda3/envs \
        -s "<comma-separated list of scaffolds>" \
        -r "<comma-separated list of reference genomes>" \
        -o <output directory>
```

#### Assess variant callers
The arguments and options of `snpeval` command:
```shell
Usage: run_benchmark.py snpeval [OPTIONS]

  SNP calling benchmark for customized dataset

Options:
  -o, --outpath PATH       The directory where to put the results and figures.
                           The path can be specified either in the CLI as
                           argument or in the config file.
  -c, --conda_prefix PATH  The prefix of conda ENV. [default: in the working
                           directory].
  -t, --threads INTEGER    The number of threads to use.  [default: 2]
  -d, --dryrun             Print the details without run the pipeline.
                           [default: False]
  -v, --vcfs TEXT          Comma-separated list of VCF files. Please quote the
                           whole parameter if there is any white space the
                           file names. The files can be specified either in
                           the CLI as argument or in the config file.
  -r, --refs TEXT          Comma-separated list of reference genome files.
                           Please quote the whole parameter if there is any
                           white space the file names. The files can be
                           specified either in the CLI as argument or in the
                           config file.
  --help                   Show this message and exit.
```

- Run the benchmarking
```shell
python3 run_benchmark.py snpeval -t 10 -c ~/miniconda3/envs \
        -v "<comma-separated list of VCF files>" \
        -r "<comma-separated list of reference genomes>" \
        -o <output directory>
```
