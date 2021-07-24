#!/usr/bin/env python3
# -- coding: utf-8 --
'''
File: run_benchmark.py
Created Date: May 16th 2019
Author: ZL Deng <dawnmsg(at)gmail.com>
---------------------------------------
Last Modified: 29th July 2019 12:03:45 pm
'''

import os
import sys
import click
import functools
import snakemake

wd = os.path.dirname(os.path.realpath(__file__))
VERSION = '0.3'


class SpecialHelpOrder(click.Group):

    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super(SpecialHelpOrder, self).__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super(SpecialHelpOrder, self).get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super(SpecialHelpOrder, self).list_commands(ctx)
        return (c[1] for c in sorted(
            (self.help_priorities.get(command, 1), command)
            for command in commands))

    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop('help_priority', 1)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(SpecialHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo("Version {}".format(VERSION))
    ctx.exit()


@click.group(cls=SpecialHelpOrder)
@click.option("--version", is_flag=True, callback=print_version,
              expose_value=False, is_eager=True, help="Print the version.")
def cli():
    pass


def common_options(f):
    options = [click.option("-d",
                            "--dryrun",
                            is_flag=True,
                            default=False,
                            show_default=True,
                            help="Print the details without run the pipeline."
                            ),
               click.option("-t",
                            "--threads",
                            type=int,
                            default=2,
                            show_default=True,
                            help="The number of threads to use."),

               click.option("-c",
                            "--conda_prefix",
                            type=click.Path(exists=True),
                            default=None,
                            help="The prefix of conda ENV. [default: in the working directory]."),
               click.option("-o",
                            "--outpath",
                            type=click.Path(),
                            default=None,
                            help="The directory where to put the results and figures. \
The path can be specified either in the CLI as argument or in the config file. \
    [default: outpath defined in the config file]")
               ]

    return functools.reduce(lambda x, opt: opt(x), options, f)


@cli.command(help_priority=1, help="Benchmarking for HCMV dataset")
@common_options
@click.option("-e",
              "--evaluation",
              required=True,
              type=click.Choice(["all", "variantcall", "assembly"]),
              help="The evaluation to run.")
@click.option("-s",
              "--slow",
              is_flag=True,
              default=False,
              show_default=True,
              help="Run the evaluation based on reads, which is very slow. \
By default, the evaluation will be based on the VCF and contig \
files provided within this software. If this parameter is on, \
this software will run all the analyses to generate outputs \
based on reads for benchmarking which is very time consuming.")
def hcmv(evaluation, dryrun=False, conda_prefix=None, slow=False, **kwargs):
    variantcall_smk = os.path.join(wd, "eval_variantcall.smk")
    assembly_smk = os.path.join(wd, "eval_assembly.smk")
    snake_kwargs = dict(runOnReads=slow)
    for arg, val in kwargs.items():
        if val != None:
            snake_kwargs[arg] = val
    if evaluation == "variantcall":
        snakes = [variantcall_smk]
    elif evaluation == "assembly":
        snakes = [assembly_smk]
    else:
        snakes = [variantcall_smk, assembly_smk]
    for snake in snakes:
        run_snake(snake, dryrun, conda_prefix, **snake_kwargs)


@cli.command(help_priority=2, help="Variants benchmark for customized dataset")
@common_options
@click.option("-v",
              "--vcfs",
              type=str,
              help="Comma-separated list of VCF files. Please quote the whole \
parameter if there is any white space the file names. The files can be \
specified either in the CLI as argument or in the config file.")
@click.option("-r",
              "--refs",
              type=str,
              help="Comma-separated list of reference genome files. Please \
quote the whole parameter if there is any white space the file names. \
(The files can be specified either in the CLI as argument or in the config file.)")
def vareval(dryrun=False, conda_prefix=None, **kwargs):
    variantcall_smk = os.path.join(wd, "eval_variant_custom.smk")
    snake_kwargs = {}
    for arg, val in kwargs.items():
        if val != None:
            snake_kwargs[arg] = val
    run_snake(variantcall_smk, dryrun, conda_prefix, **snake_kwargs)


@cli.command(help_priority=3, help="Assembly benchmark for customized dataset")
@common_options
@click.option("-s",
              "--scaffolds",
              type=str,
              help="Comma-separated list of scaffold files. Please quote the \
whole parameter if there is any white space the file names. \
The files can be specified either in the CLI as argument or in the config file.")
@click.option("-r",
              "--refs",
              type=str,
              help="Comma-separated list of reference genome files. Please \
quote the whole parameter if there is any white space in the file names. \
(The files can be specified either in the CLI as argument or in the config file.)")
def asmeval(dryrun=False, threads=2, conda_prefix=None, **kwargs):
    #snpcall_smk = os.path.join(wd, "evaluate_snpcall_customize.smk")
    assembly_smk = os.path.join(wd, "eval_assembly_custom.smk")
    snake_kwargs = {}
    for arg, val in kwargs.items():
        if val != None:
            snake_kwargs[arg] = val
    run_snake(assembly_smk, dryrun, conda_prefix, **snake_kwargs)


def run_snake(snake, dryrun=False, conda_prefix=None, **kwargs):
    try:
        # Unlock the working directory
#         unlocked = snakemake.snakemake(
#             snakefile=snake,
#             # unlock=False,
#             unlock=True,
#             workdir=wd,
#             config=kwargs
#         )
#         if not unlocked:
#             raise Exception('Could not unlock the working directory!')

        # Start the snakemake pipeline
        success = snakemake.snakemake(
            snakefile=snake,
            restart_times=1,
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
