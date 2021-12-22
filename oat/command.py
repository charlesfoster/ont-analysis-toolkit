#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last edited on Thur December 02 2021

@author: Dr Charles Foster
"""

#==============================================================#
# %% Import Modules                                               #
#==============================================================#
import os
import shutil
import sys
import psutil
from datetime import date
import argparse
from argparse import RawTextHelpFormatter
from oat.scripts.helper_functions import find_runDirs, check_arguments, initiate_colorlog, printc
import pandas as pd
from psutil import virtual_memory
from oat import __version__
#=============================================================#
# %%GLOBAL VARIABLES AND DEPENDENCIES                            #
#=============================================================#
global main_dir
TODAY = date.today().strftime("%Y-%m-%d")
hash_line = "#"*46

thisdir = os.path.abspath(os.path.dirname(__file__))
main_dir = os.path.split(thisdir)[0]
snakefile = os.path.join(thisdir, "scripts", "analysis_module.smk")

# ============================================================#
# Define functions                                            #
# ============================================================#

def bytesto(bytes, to, bsize=1024):
    """convert bytes to megabytes, etc.
    sample code:
        print('mb= ' + str(bytesto(314575262000000, 'm')))
    sample output:
        mb= 300002347.946
    """

    a = {"k": 1, "m": 2, "g": 3, "t": 4, "p": 5, "e": 6}
    r = float(bytes)
    for i in range(a[to]):
        r = r / bsize

    return r

#def retrieve_globals

# %% main
def main(sysargs=sys.argv[1:]):
    print(
        """\n\033[95m 
                                        ,d
                                        88
               ,adPPYba,   ,adPPYYba, MM88MMM
              a8"     "8a ""     `Y8    88   
              8b       d8 ,adPPPPP88    88   
              "8a,   ,a8" 88,    ,88    88,  
               `"YbbdP"'  `"8bbdP"Y8   "Y888
        
        OAT: ONT Analysis Toolkit (version {})\033[0m
    """.format(
            __version__
        )
    )
    
    max_mem = round(bytesto(virtual_memory().available, "m"))
    
    parser = argparse.ArgumentParser(
        description="A pipeline for sequencing and analysis of viral genomes using an ONT MinION",
        usage="""oat [options] <samples_file> """,
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "samples_file", nargs="*", help="Path to file with sample metadata (.csv format). See example spreadsheet for minimum necessary information."
    )
    parser.add_argument(
        "-b",
        "--barcode_kit",
        action="store",
        required=False,
        default="SQK-RBK004",
        help="""
            Barcode kit that you used: necessary for demultiplexing.
            
            Common options include:
            - rapid 12-barcode kit (SQK-RBK004)
            - rapid 96-barcode kit (SQK-RBK110-96)
            - native 12-barcode kit (EXP-NBD104)
            - native 12-barcode expansion kit 13-24 (EXP-NBD114)
            - native 96-barcode kit (EXP-NBD196)
            
            Please ensure your version of guppy_barcoder supports the kit name. Try: 'guppy_barcoder --print_kits'
            
            Default: {}
            """.format(
            'SQK-RBK004'
        ),
    )
    parser.add_argument(
        "-c",
        "--consensus_freq",
        action="store",
        required=False,
        default=float(0.8),
        help="""
            Variant allele frequency threshold for a variant to be incorporated into consensus genome.
            Variants below this frequency will be incorporated with an IUPAC ambiguity. 
            Default: {}
            """.format(
            float(0.80)
        ),
        metavar="<float>",
    )
    parser.add_argument(
        "-d",
        "--demultiplexed",
        action="store_true",
        default=False,
        required=False,
        help="""
            Reads already demultiplexed using guppy_barcoder into '/var/lib/minknow/data/<run_name>'.
            By default, assumes reads need to be demultiplexed and reads are demultiplexed into the output directory.
            """,
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        required=False,
        help="Force overwriting of completed files in snakemake analysis (Default: files not overwritten)",
    )
    parser.add_argument(
        "-n", "--dry_run", action="store_true", default=False, help="Dry run only"
    )
    parser.add_argument(
        "-m", 
        "--module", 
        action="store", 
        default="all", 
        required=False,
        help="Pipeline module to run: 'rampart', 'analysis' or 'all' (rampart followed by analysis). Default: 'all'", metavar='rampart | analysis | all'
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Output directory. Default: {} + 'run_name' from samples spreadsheet".format(
            os.path.join(os.getcwd(), "analysis_results")
        ),
    )
    parser.add_argument(
        "--rampart_outdir",
        action="store",
        help="Output directory. Default: {}".format(
            os.path.join(os.getcwd(), "rampart_files")
        ),
    )
    parser.add_argument(
        "-p",
        "--print_dag",
        action="store_true",
        default=False,
        help="Save directed acyclic graph (DAG) of workflow to outdir",
    )
    parser.add_argument(
        "-r",
        "--reference",
        action="store",
        required=False,
        default="MN908947.3",
        help="""
            Reference genome to use. Supported references:
                - 'MN908947.3' (SARS-CoV-2)
                - 'NC_006273.2' (CMV Merlin) (planned support, but not yet implemented)
            Other references can be used, but the corresponding assembly (fasta) and annotation (gff3 from Ensembl) must be added. See README.
            Default: {0}
            """.format(
            "MN908947.3"
        ),
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        help="Number of threads to use",
        default=psutil.cpu_count(logical=True),
        metavar="<int>",
    )
    parser.add_argument(
        "-v",
        "--variant_caller",
        action="store",
        help="Variant caller to use. Choices: 'medaka-longshot'. Default: 'medaka-longshot'",
        default="medaka-longshot",
    )
    parser.add_argument(
        "--create_envs_only",
        action="store_true",
        default=False,
        help="Create conda environments for snakemake analysis, but do no further analysis. Useful for initial pipeline setup. Default: False",
    )
    parser.add_argument(
        "--snv_min",
        action="store",
        help="Minimum variant allele frequency for an SNV to be kept Default: 0.2",
        default=float(0.4),
    )
    parser.add_argument(
        "--delete_reads", action="store_true", help="Delete demultiplexed reads after analysis", default=False
    )
    parser.add_argument(
        "--redo_analysis", action="store_true", help="Delete entire analysis output directory and contents for a fresh run", default=False
    )
    parser.add_argument(
        "--version",
        action="version",
        version="covid illumina pipeline snakemake edition:  0.1.0",
        #version=f"covid illumina pipeline snakemake edition:  {__version__}",
    )
    parser.add_argument(
        "--minknow_data",
        action="store",
        default="/var/lib/minknow/data",
        help="Location of MinKNOW data root. Default: {}".format(
            "/var/lib/minknow/data"
            ),
    )
    parser.add_argument(
        "--max_memory",
        action="store",
        help="Maximum memory (in MB) that you would like to provide to snakemake. Default: {}MB".format(
            max_mem
        ),
        metavar="<int>",
    )
    parser.add_argument("--quiet", action="store_true", help="Stop printing of snakemake commands to screen.")
    parser.add_argument("--report", action="store_true", help="Generate report (currently non-functional).")
    
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    args = parser.parse_args()

    ### end parsing of command args ###

    if args.print_dag:      
        args.module = "analysis"
    args.module = args.module.upper()
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    
    minknow_dir = args.minknow_data
    if not os.path.exists(minknow_dir):
        print(
            "#####\n\033[91mError\033[0m: The 'minknow_data' path appears to be incorrect. Check.\n#####\n"
        )
        sys.exit(1)

    variable_dict = vars(args)
    variable_dict["max_mem"] = max_mem
    initiate_colorlog(variable_dict, main_dir)
    check_arguments(variable_dict, args)
    find_runDirs(variable_dict, main_dir, minknow_dir)
    my_log = variable_dict["my_log"]
    
    if args.module == 'RAMPART':
        from oat.scripts.rampart_module import rampart_json, rampart_run
        if args.rampart_outdir:
            rampart_outdir = ''.join(args.rampart_outdir)
        else:   
            rampart_outdir = os.path.join(os.getcwd(), "rampart_files")
        variable_dict['rampart_outdir'] = rampart_outdir
        rampart_json(variable_dict)
        rampart_run(variable_dict)
        shutil.move(variable_dict['logfile'], os.path.join(main_dir,TODAY+'_'+variable_dict['run_name']+'_RAMPART.log'))
        printc("\n Pipeline complete\n", "HEADER")
    elif args.module == 'ANALYSIS' or args.module == 'ALL':
        if args.outdir:
            analysis_outdir = ''.join(args.outdir)
        else:   
            analysis_outdir = os.path.join(os.getcwd(), "analysis_results",variable_dict["run_name"])
        variable_dict['outdir'] = analysis_outdir
        if args.module == 'ALL':
            #first run rampart
            from oat.scripts.rampart_module import rampart_json, rampart_run
            if args.rampart_outdir:
                rampart_outdir = ''.join(args.rampart_outdir)
            else:   
                rampart_outdir = os.path.join(os.getcwd(), "rampart_files")
            variable_dict['rampart_outdir'] = rampart_outdir
            rampart_json(variable_dict)
            rampart_run(variable_dict)
            final_log_name = os.path.join(main_dir,TODAY+'_'+variable_dict['run_name']+'_ALL.log')
        else:
            final_log_name = os.path.join(main_dir,TODAY+'_'+variable_dict['run_name']+'_ANALYSIS.log')
        if args.redo_analysis:
            my_log.info("You have chosen to redo the analyses for {0}".format(variable_dict["run_name"]))
            decision = input("Are you sure you want to redo the analysis? This option will delete all analysis data for {0}. Type 'yes' to continue, or anything else to abort.\n> ".format(variable_dict["run_name"]))
            if decision.upper() == "YES":
                decision = input("Is this the correct folder to delete? {0}. Type 'yes' to continue, or anything else to abort.\n> ".format(variable_dict["outdir"]))
                if decision.upper() == "YES":
                    shutil.rmtree(variable_dict["outdir"])
                else:
                    sys.exit("You have chosen not to delete the old data after all. Quitting.")
            else:
                sys.exit("You have chosen not to redo the analysis after all. Quitting.")
        if not os.path.exists(variable_dict['outdir']):
           os.makedirs(variable_dict['outdir'])
       # prepare parameters for analysis
        length_params= pd.read_csv(os.path.join(thisdir,"protocols","length_params.csv")).set_index("name")
        variable_dict['min_len'] = length_params.loc[variable_dict['protocol'], 'min_length']
        variable_dict['max_len'] = length_params.loc[variable_dict['protocol'], 'max_length']     
        variable_dict["run_data"].to_csv(os.path.join(variable_dict["outdir"],"metadata.csv"))
        
        if not args.demultiplexed:
            from oat.scripts.demux_and_filter import demultiplex_reads, filter_reads
            demultiplex_reads(variable_dict)
            filter_reads(variable_dict)
        else:
            from oat.scripts.demux_and_filter import relocate_and_filter_reads
            relocate_and_filter_reads(variable_dict)        
        
        if args.print_dag:
            flat_config = []
            variable_dict["run_data"] = ''.join(args.samples_file)
            del variable_dict["my_log"]
            del variable_dict["barcodes_used"]
            for key in variable_dict:
                flat_config.append(key + "=" + str(variable_dict[key]))
            flat_config = " ".join(flat_config)
            cmd = 'snakemake -j1 -s {0} --quiet --config {1} --dag | grep -v "No negative control samples detected" | dot -Tpdf > {2}'.format(
                snakefile, 
                flat_config,
                os.path.join(variable_dict["outdir"], "workflow_DAG.pdf")
            )
            os.system(cmd)
            my_log.info("Saving directed acyclic graph of to {}".format(os.path.join(variable_dict["outdir"], "workflow_DAG.pdf")))
            my_log.info("To run the analysis, use the same command but omit '-p' / '--print_dag'")
            printc("\n Pipeline complete\n", "HEADER")
            sys.exit(0)

        # time for some snakemake action
        import snakemake
        my_log.info("Running analysis pipeline using snakemake")

        #check if user only wants to create conda environments     
        if args.create_envs_only:
            status = snakemake.snakemake(
                snakefile,
                use_conda=True,
                conda_frontend="mamba",
                conda_create_envs_only=True,
                dryrun=False,
                printshellcmds=True,
                forceall=False,
                force_incomplete=True,
                config=variable_dict,
                quiet=False,
                cores=args.threads,
                lock=False,
            )
            if status:
                print(
                    "\033[Environments created!\033[0m\n"
                    )
                sys.exit(0)
            else:
                my_log.error("Something went wrong with conda environment creation. Investigate.")
                sys.exit(1)
        
        # if user does not want that option
        elif variable_dict["quiet"]:
            status = snakemake.snakemake(
                snakefile,
                use_conda=True,
                conda_frontend="mamba",
                dryrun=args.dry_run,
                printshellcmds=False,
                forceall=args.force,
                force_incomplete=True,
               # resources=max_mem,
                config=variable_dict,
                quiet=True,
                cores=args.threads,
                lock=False,
            )
        else:
            print("\n**** CONFIG ****")
            for k in variable_dict:
                if k not in ["my_log", "run_data", "password"]:
                    print(k + ": ", variable_dict[k])
            print("")
            status = snakemake.snakemake(
                snakefile,
                use_conda=True,
                conda_frontend="mamba",
                dryrun=args.dry_run,
                printshellcmds=True,
                forceall=args.force,
                force_incomplete=True,
                config=variable_dict,
                quiet=False,
                cores=args.threads,
                lock=False,
            )
        shutil.move(variable_dict['logfile'], final_log_name)
        
        my_log.info("Analysis module complete")
        my_log.info("Final results summary: {}".format(os.path.join(variable_dict['outdir'],variable_dict['run_name']+'_qc.csv')))
        printc("\n Pipeline complete\n", "HEADER")
##########
# %% run analysis                                                     #

if __name__ == '__main__':
    main()
    
