#!/usr/bin/env python3
import os
import sys
import argparse
from argparse import RawTextHelpFormatter

wd = os.path.dirname(os.path.realpath(__file__))

def start_evaluation(snake, dryrun=False, threads=2, conda_prefix=None):
    try:
        import snakemake
        
        ## Unlock the working directory
        unlocked = snakemake.snakemake(
            snakefile=snake,
            #unlock=False,
            unlock=True,
            workdir=wd
        )
        if not unlocked:
            raise Exception('Could not unlock the working directory!')

        ## Start the snakemake pipeline
        success=snakemake.snakemake(
            snakefile=snake,
            restart_times=3, 
            cores=threads,
            workdir=wd,
            use_conda=True,
            conda_prefix=conda_prefix,
            dryrun= dryrun,
            summary=True,
            printshellcmds= True,
            force_incomplete= True
            )
        if not success:
            raise Exception('Snakemake pipeline failed!')
    except Exception as e:
        from datetime import datetime
        print('ERROR ({})'.format(snake))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ',timespec= 'minutes'),
            e))
        raise RuntimeError(e)
    except :
        from datetime import datetime
        print('ERROR ({})'.format(snake))
        print('{}\t{}\n'.format(
            datetime.now().isoformat(' ',timespec= 'minutes'),
            sys.exc_info()))
        raise RuntimeError('Unknown problem occured when lauching Snakemake!')

if __name__ == "__main__":
    usage = r'''
    run_benchmark.py --- run the benchmarking for assembly and SNPs calling


    Usage:
    python run_benchmark.py [-d] [-t threads] <all|snpcall|assembly>

    '''

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=RawTextHelpFormatter,
                                     epilog="   by ZL Deng")

    parser.add_argument("evaluation",
            type=str,
            choices=["all", "snpcall", "assembly"],
                        help="the evaluation to run")
    parser.add_argument("-d",
                    "--dryrun",
                    dest="dryrun",
                    default=False, 
            action="store_true",
                    help="Print the details without run the pipeline")

    parser.add_argument("-t", 
            "--threads", 
            type=int, 
            default=2,
            help="The number of threads to use, default: 2")
    parser.add_argument("-c", 
            "--conda_prefix", 
            type=str, 
            default=None,
            help="The prefix of conda ENV")
    args = parser.parse_args()
    snpcall_smk = os.path.join(wd, "evaluate_snpcall.smk")
    assembly_smk = os.path.join(wd, "evaluate_assembly.smk")
    if args.evaluation == "snpcall":
        start_evaluation(snpcall_smk, args.dryrun, args.threads, args.conda_prefix)
    elif args.evaluation == "assembly":
        start_evaluation(assembly_smk, args.dryrun, args.threads, args.conda_prefix)
    else:
        start_evaluation(snpcall_smk, args.dryrun, args.threads, args.conda_prefix)
        start_evaluation(assembly_smk, args.dryrun, args.threads, args.conda_prefix)

