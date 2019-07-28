#!/usr/bin/env python3
# -- coding: utf-8 --
import os
import sys
import click
import functools
import snakemake

wd = os.path.dirname(os.path.realpath(__file__))


@click.group()
def cli():
    pass


def common_options(f):
    options = [click.option("-d",
                            "--dryrun",
                            is_flag=True,
                            default=False,
                            help="Print the details without run the pipeline"
                            ),
               click.option("-t",
                            "--threads",
                            type=int,
                            default=2,
                            help="The number of threads to use, default: 2"),

               click.option("-c",
                            "--conda_prefix",
                            type=click.Path(exists=True),
                            default=None,
                            help="The prefix of conda ENV [default: in the working directory]"),
               click.option("-o",
                            "--outpath",
                            type=click.Path(),
                            default=None,
                            help="The directory where to put the results and figures")
               ]

    return functools.reduce(lambda x, opt: opt(x), options, f)


@cli.command()
@common_options
@click.option("-e",
              "--evaluation",
              required=True,
              type=click.Choice(["all", "snpcall", "assembly"]),
              help="The evaluation to run.")
@click.option("-s",
              "--slow",
              is_flag=True,
              default=False,
              help="Run the evaluation based on reads, which is very slow. \
By default, the evaluation will be based on the VCF and contig \
files provided within this software. If this parameter is on, \
this software will run all the analyses to generate outputs \
based on reads for benchmarking which is very time consuming.")
def hcmv(evaluation, dryrun=False, conda_prefix=None, slow=False, **kwargs):
    snpcall_smk = os.path.join(wd, "evaluate_snpcall.smk")
    assembly_smk = os.path.join(wd, "evaluate_assembly.smk")
    snake_kwargs = dict(runOnReads=slow)
    for arg, val in kwargs.items():
        if val != None:
            snake_kwargs[arg] = val
    if evaluation == "snpcall":
        snakes = [snpcall_smk]
    elif evaluation == "assembly":
        snakes = [assembly_smk]
    else:
        snakes = [snpcall_smk, assembly_smk]
    for snake in snakes:
        run_snake(snake, dryrun, conda_prefix, **snake_kwargs)


@cli.command()
@common_options
@click.option("-v",
              "--vcfs",
              type=str,
              help="Comma-separated list of VCF files. Please quote the whole \
parameter if there is any white space the file names")
@click.option("-r",
              "--refs",
              type=str,
              help="Comma-separated list of reference genome files. Please \
quote the whole parameter if there is any white space the file names")
def snpeval(dryrun=False, conda_prefix=None, **kwargs):
    snpcall_smk = os.path.join(wd, "evaluate_snpcall_customize.smk")
    snake_kwargs = {}
    for arg, val in kwargs.items():
        if val != None:
            snake_kwargs[arg] = val
    run_snake(snpcall_smk, dryrun, conda_prefix, **snake_kwargs)


@cli.command()
@common_options
@click.option("-s",
              "--scaffolds",
              type=str,
              help="Comma-separated list of scaffold files. Please quote the \
whole parameter if there is any white space the file names")
@click.option("-r",
              "--refs",
              type=str,
              help="Comma-separated list of reference genome files. Please \
quote the whole parameter if there is any white space the file names")
def asmeval(dryrun=False, threads=2, conda_prefix=None, **kwargs):
    #snpcall_smk = os.path.join(wd, "evaluate_snpcall_customize.smk")
    assembly_smk = os.path.join(wd, "evaluate_assembly_customize.smk")
    snake_kwargs = {}
    for arg, val in kwargs.items():
        if val != None:
            snake_kwargs[arg] = val
    run_snake(assembly_smk, dryrun, conda_prefix, **snake_kwargs)


def run_snake(snake, dryrun=False, conda_prefix=None, **kwargs):
    try:
        # Unlock the working directory
        unlocked = snakemake.snakemake(
            snakefile=snake,
            # unlock=False,
            unlock=True,
            workdir=wd,
            config=kwargs
        )
        if not unlocked:
            raise Exception('Could not unlock the working directory!')

        # Start the snakemake pipeline
        success = snakemake.snakemake(
            snakefile=snake,
            restart_times=3,
            cores=kwargs.get("threads", 2),
            workdir=wd,
            use_conda=True,
            conda_prefix=conda_prefix,
            dryrun=dryrun,
            printshellcmds=True,
            force_incomplete=True,
            config=kwargs,
        )
        if not success:
            raise Exception('Snakemake pipeline failed!')
    except Exception as e:
        from datetime import datetime
        print('ERROR ({})'.format(snake))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ', timespec='minutes'),
            e))
        raise RuntimeError(e)
    except:
        from datetime import datetime
        print('ERROR ({})'.format(snake))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ', timespec='minutes'),
            sys.exc_info()))
        raise RuntimeError(
            'Unknown problem occured when lauching Snakemake!')


if __name__ == "__main__":
    cli()
