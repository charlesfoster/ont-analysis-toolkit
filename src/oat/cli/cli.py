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
import re
import shutil
import sys
import psutil
from datetime import date
import argparse
from argparse import RawTextHelpFormatter
from oat.cli.scripts.helper_functions import find_runDirs, check_arguments, initiate_colorlog, printc, check_prior_lineages
import pandas as pd
from psutil import virtual_memory
from oat import __version__
import GPUtil
import textwrap as _textwrap

#=============================================================#
# %%GLOBAL VARIABLES AND DEPENDENCIES                            #
#=============================================================#
global main_dir
TODAY = date.today().strftime("%Y-%m-%d")
hash_line = "#"*46

thisdir = os.path.abspath(os.path.dirname(__file__))
main_dir = os.path.split(thisdir)[0]
snakefile = os.path.join(thisdir, "scripts", "analysis_module.smk")
alternate_snakefile = os.path.join(thisdir, "scripts", "alternate_analysis.smk")


# Define the wrapping of help text: https://stackoverflow.com/questions/35917547/python-argparse-rawtexthelpformatter-with-line-wrap
os.environ['COLUMNS'] = "120"
class PreserveWhiteSpaceWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def __add_whitespace(self, idx, iWSpace, text):
        if idx == 0:
            return text
        return (" " * iWSpace) + text

    def _split_lines(self, text, width):
        textRows = text.splitlines()
        for idx,line in enumerate(textRows):
            search = re.search('\s*[0-9\-]{0,}\.?\s*', line)
            if line.strip() == "":
                textRows[idx] = " "
            elif search:
                lWSpace = search.end()
                lines = [self.__add_whitespace(i,lWSpace,x) for i,x in enumerate(_textwrap.wrap(line, width))]
                textRows[idx] = lines

        return [item for sublist in textRows for item in sublist]

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
    # initial resource wrangling
    max_mem = round(bytesto(virtual_memory().available, "m"))
    num_gpu = len(GPUtil.getGPUs())

    parser = argparse.ArgumentParser(
        description="A pipeline for sequencing and analysis of viral genomes using an ONT MinION",
        usage="""oat [options] <samples_file> """,
        formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
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
            - new rapid 24-barcode kit (SQK-RBK114-24)
            - new rapid 96-barcode kit (SQK-RBK114-96)
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
        default=float(0),
        help="""
            Variant allele frequency threshold for a variant to be incorporated into consensus genome.
            Variants below this frequency will be incorporated with an IUPAC ambiguity.
            Set to 0 to incorporate the majority or most common base.
            Note: currently do not recommend anything except the default - debugging.
            Default: {}
            """.format(
            float(0)
        ),
        metavar="<float>",
    )
    parser.add_argument(
        "-i",
        "--indel_freq",
        action="store",
        required=False,
        default=float(0.40),
        help="""
            Variant allele frequency threshold for an indel variant to be incorporated into consensus genome.
            Variants below this frequency will not be incorporated.
            Set to 0 to incorporate the majority or most common base.
            Default: {}
            """.format(
            float(0.40)
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
        help="Output directory. Default: {} + 'organism_name' + 'run_name' from samples spreadsheet".format(
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
        help="Variant caller to use. Choices: 'clair3','medaka'. Default: 'clair3'",
        default="clair3",
    )
    parser.add_argument(
        "--create_envs_only",
        action="store_true",
        default=False,
        help="Create conda environments for snakemake analysis, but do no further analysis. Useful for initial pipeline setup. Default: False",
    )
    parser.add_argument(
        "-s",
        "--snv_min_freq",
        action="store",
        help="Minimum variant allele frequency for an SNV to be kept. Default: 0.2",
        default=float(0.2),
    )
    parser.add_argument(
        "--guppy_model",
        action="store",
        help="Model used within guppy for basecalling - needed for medaka analyses. 'GUPPY_MODEL' environmental model is used, if set, otherwise default: r941_min_high_g360",
        default="r941_min_high_g360",
    )
    parser.add_argument(
        "--clair3_model",
        action="store",
        help="Path to where clair3 model is located. 'CLAIR3_MODEL' environmental model is used, if set, otherwise singularity default: /opt/models/r941_prom_hac_g360+g422",
        default="/opt/models/r941_prom_hac_g360+g422",
    )
    parser.add_argument(
        "-M",
        "--min_depth",
        action="store",
        help="Minimum depth for (1) an SNV to be kept; and (2) consensus genome generation. Default: 20",
        default=int(20),
    )
    parser.add_argument(
        "--alternate_analysis",
        action="store_true",
        help="Run an alternate analysis to generate consensus genomes based on different input parameters. Use after initial run.",
        default=False,
    )
    parser.add_argument(
        "--alt_cov_max",
        action="store",
        help="Samples with a genome coverage alt_cov_min <= x < alt_cov_max will be chosen for an alternate analysis. ONLY USED WITH ALTERNATE ANALYSIS OPTION. Default: 80",
        default=float(80),
    )
    parser.add_argument(
        "--alt_cov_min",
        action="store",
        help="Samples with a genome coverage alt_cov_min <= x < alt_cov_max will be chosen for an alternate analysis. ONLY USED WITH ALTERNATE ANALYSIS OPTION. Default: 40",
        default=float(40),
    )
    parser.add_argument(
        "--delete_reads", action="store_true", help="Delete demultiplexed reads after analysis", default=False
    )
    parser.add_argument(
        "--redo_analysis", action="store_true", help="Delete entire analysis output directory and contents for a fresh run", default=False
    )
    parser.add_argument(
        "--additional_nanoq", action="store", help="Additional results for nanoq besides default length filtering", default=""
    )
    parser.add_argument(
        "--skip_clipping", action="store_true", help="Skip clipping of amplicon primers", default=False
    )
    parser.add_argument(
        "--no_barcodes", action="store_true", help="No barcodes were used during library prep (or pretend this is so)", default=False
    )
    parser.add_argument(
        "--version",
        action="version",
        #version="covid illumina pipeline snakemake edition:  0.1.0",
        version=f"OAT: ONT Analysis Toolkit:  v{__version__}",
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
        "--no_update",
        action="store_true",
        default=False,
        help="Disable updating of container versions for SARS-COV-2 analysis. Default: {}".format(
            False
            ),
    )
    parser.add_argument(
        "--list_protocols",
        action="store_true",
        default=False,
        help="List available protocols and exit",
    )
    parser.add_argument(
        "-MM",
        "--max_memory",
        action="store",
        help="Maximum memory (in MB) that you would like to provide to snakemake. Default: {}MB".format(
            max_mem
        ),
        metavar="<int>",
    )
    parser.add_argument("--quiet", action="store_true", help="Stop printing of Snakemake commands to screen.")
    parser.add_argument("--report", action="store_true", help="Generate report (currently minimally functional).")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    args = parser.parse_args()

    ### end parsing of command args ###
    args.module = args.module.upper()
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)

    ### check if listing protocols
    if args.list_protocols:
        from oat.cli.scripts.helper_functions import list_protocols
        list_protocols()
        sys.exit(0)

    minknow_dir = args.minknow_data
    if not os.path.exists(minknow_dir):
        print(
            "#####\n\033[91mError\033[0m: The 'minknow_data' path appears to be incorrect. Check.\n#####\n"
        )
        sys.exit(1)

    variable_dict = vars(args)

    if args.max_memory:
        max_mem = int(args.max_memory)
    resources_dict = {'mem_mb':max_mem,
                     'gpu':num_gpu}

    variable_dict["resources"] = resources_dict

    initiate_colorlog(variable_dict, os.getcwd())
    check_arguments(variable_dict, args)
    find_runDirs(variable_dict, main_dir, minknow_dir)
    my_log = variable_dict["my_log"]


    if args.module == 'RAMPART':
        from oat.cli.scripts.rampart_module import rampart_json, rampart_run, rampart_watchdog
        if args.rampart_outdir:
            rampart_outdir = ''.join(args.rampart_outdir)
        else:
            rampart_outdir = os.path.join(os.getcwd(), "rampart_files")
        variable_dict['rampart_outdir'] = rampart_outdir
        rampart_json(variable_dict)
        if variable_dict['no_barcodes']:
            rampart_watchdog(variable_dict)
        else:
            rampart_run(variable_dict)
        #set up correct log destination
        if args.outdir:
            logdir = args.outdir
        else:
            logdir = os.getcwd()
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        shutil.move(variable_dict['logfile'], os.path.join(logdir,TODAY+'_'+variable_dict['run_name']+'_RAMPART.log'))
        printc("\n Pipeline complete\n", "HEADER")
    elif args.module == 'ANALYSIS' or args.module == 'ALL' or args.module == 'ALTERNATE':
        if args.module == 'ALL':
            #first run rampart
            from oat.cli.scripts.rampart_module import rampart_json, rampart_run
            if args.rampart_outdir:
                rampart_outdir = ''.join(args.rampart_outdir)
            else:
                rampart_outdir = os.path.join(os.getcwd(), "rampart_files")
            variable_dict['rampart_outdir'] = rampart_outdir
            rampart_json(variable_dict)
            rampart_run(variable_dict)
            final_log_name = os.path.join(variable_dict['outdir'],TODAY+'_'+variable_dict['run_name']+'_ALL.log')
        elif args.module =="ANALYSIS":
            final_log_name = os.path.join(variable_dict['outdir'],TODAY+'_'+variable_dict['run_name']+'_ANALYSIS.log')
        if args.redo_analysis and not args.alternate_analysis:
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
        elif args.redo_analysis and args.alternate_analysis:
            sys.exit("You cannot redo the analysis from scratch and run the 'alternate analysis' at the same time. Quitting.")

        if not os.path.exists(variable_dict['outdir']):
           os.makedirs(variable_dict['outdir'])

        # prepare parameters for analysis
        length_params= pd.read_csv(os.path.join(thisdir,"protocols","length_params.csv")).set_index("name")
        variable_dict['min_len'] = length_params.loc[variable_dict['protocol'], 'min_length']
        variable_dict['max_len'] = length_params.loc[variable_dict['protocol'], 'max_length']
        if args.alternate_analysis:
            variable_dict["run_data"].to_csv(os.path.join(variable_dict["outdir"],"metadata.alternate.csv"), index=False)
        else:
            variable_dict["run_data"].to_csv(os.path.join(variable_dict["outdir"],"metadata.csv"), index=False)

        # run the alternate analysis
        if args.alternate_analysis:
            final_log_name = os.path.join(variable_dict['outdir'],TODAY+'_'+variable_dict['run_name']+'_ALTERNATE_ANALYSIS.log')
            # time for some snakemake action
            import snakemake
            my_log.info("Running alternate analysis pipeline using snakemake")
            my_log.info("If you have not previously run the pipeline 'normally', this analysis will fail.")

            # find previous qc file and read it
            previous_qc = os.path.join(variable_dict['outdir'],variable_dict['run_name']+'_qc.csv')
            if not os.path.exists(previous_qc):
                sys.exit("Could not locate QC file for initial analysis. Investigate.")
            else:
                df = pd.read_csv(previous_qc)

            # select samples with low coverage
            if 'coverage' in df.columns:
                low_cov_samples = list(df[df['coverage'].between(float(args.alt_cov_min), float(args.alt_cov_max))]['id'])
            else:
                # deal with legacy colname
                cov_col = ''.join(list(df.filter(regex='ref_cov_*').columns))
                df['coverage'] = [float(x) for x in list(df[cov_col])]
                low_cov_samples = list(df[df['coverage'].between(float(args.alt_cov_min), float(args.alt_cov_max))]['id'])

            # die if no low cov samples
            if len(low_cov_samples) == 0:
                sys.exit("No samples needing alternate analysis were identified based on your input parameters. Quitting.")

            variable_dict['alternate_isolates'] = low_cov_samples

            # if a SARS-CoV-2 analysis, delete previous lineage files to trigger re-run
            if os.path.basename(variable_dict["reference"]) == "MN908947.3.fasta":
                check_prior_lineages(variable_dict)

            # run alternate analysis
            print("\n**** CONFIG ****")
            for k in variable_dict:
                if k not in ["my_log", "run_data", "password"]:
                    print(k + ": ", variable_dict[k])
            print("")
            status = snakemake.snakemake(
                alternate_snakefile,
                use_conda=True,
                use_singularity=True,
                singularity_args=variable_dict['singularity_args'],
                conda_frontend="mamba",
                dryrun=args.dry_run,
                printshellcmds=True,
                forceall=args.force,
                force_incomplete=True,
                resources=variable_dict['resources'],
                list_resources=False,
                config=variable_dict,
                quiet=False,
                cores=int(args.threads),
                lock=False,
            )

            shutil.move(variable_dict['logfile'], final_log_name)
            if status:
                my_log.info("Alternate analysis complete")
                my_log.info("Final results summary: {}".format(os.path.join(variable_dict['outdir'],variable_dict['run_name']+'_qc.alternate.csv')))
                printc("\n Pipeline complete\n", "HEADER")
                sys.exit(0)
            else:
                my_log.error("Something went wrong with the alternate analysis. Investigate.")
                sys.exit(1)


        if not args.demultiplexed:
            from oat.cli.scripts.demux_and_filter import demultiplex_reads, filter_reads
            demultiplex_reads(variable_dict)
            my_log.info("Filtering reads")
            filter_reads(variable_dict)
        else:
            from oat.cli.scripts.demux_and_filter import relocate_and_filter_reads
            my_log.info("Relocating and filtering reads")
            relocate_and_filter_reads(variable_dict)

        if args.print_dag:
            flat_config = []
            variable_dict["run_data"] = ''.join(args.samples_file)
            del variable_dict["my_log"]
            del variable_dict["barcodes_used"]
            del variable_dict["resources"]
            del variable_dict["singularity_args"]
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
        elif args.report:
            my_log.info("Generating snakemake report")
            my_log.warning("Reporting is currently minimal/not very useful")
            import snakemake
            if not os.path.exists(variable_dict['outdir']):
                os.makedirs(variable_dict['outdir'])
            status = snakemake.snakemake(
                snakefile,
                report=os.path.join(variable_dict['outdir'], "pipeline_report.html"),
                use_conda=True,
                use_singularity=True,
                singularity_args=variable_dict['singularity_args'],
                conda_frontend="mamba",
                dryrun=args.dry_run,
                printshellcmds=True,
                forceall=args.force,
                force_incomplete=True,
                resources=variable_dict['resources'],
                config=variable_dict,
                quiet=True,
                cores=1,
                lock=False,
            )
            if status:
                my_log.info("View report: {}".format(os.path.join(variable_dict['outdir'], "pipeline_report.html")))
                print(
                    "\033[92m\nReport created!\033[0m\n"
                    )
                sys.exit(0)
            else:
                my_log.error("Something went wrong with report creation. Investigate.")
                sys.exit(1)

        # time for some snakemake action
        import snakemake
        my_log.info("Running analysis pipeline using snakemake")

        # if a SARS-CoV-2 analysis, delete previous lineage files to trigger re-run
        if os.path.basename(variable_dict["reference"]) == "MN908947.3.fasta":
            check_prior_lineages(variable_dict)

        #check if user only wants to create conda environments
        if args.create_envs_only:
            status = snakemake.snakemake(
                snakefile,
                use_conda=True,
                use_singularity=True,
                singularity_args=variable_dict['singularity_args'],
                conda_frontend="mamba",
                conda_create_envs_only=True,
                dryrun=False,
                printshellcmds=True,
                forceall=False,
                force_incomplete=True,
                config=variable_dict,
                quiet=False,
                cores=int(args.threads),
                lock=False,
            )
            if status:
                print(
                    "\033[92m\nEnvironments created!\033[0m\n"
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
                use_singularity=True,
                singularity_args=variable_dict['singularity_args'],
                conda_frontend="mamba",
                dryrun=args.dry_run,
                printshellcmds=False,
                forceall=args.force,
                force_incomplete=True,
                resources=variable_dict['resources'],
                config=variable_dict,
                quiet=True,
                cores=int(args.threads),
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
                use_singularity=True,
                singularity_args=variable_dict['singularity_args'],
                conda_frontend="mamba",
                dryrun=args.dry_run,
                printshellcmds=True,
                forceall=args.force,
                force_incomplete=True,
                resources=variable_dict['resources'],
                list_resources=False,
                config=variable_dict,
                quiet=False,
                cores=int(args.threads),
                lock=False,
            )
        shutil.move(variable_dict['logfile'], final_log_name)

        if status:
            my_log.info("Analysis module complete")
            my_log.info("Final results summary: {}".format(os.path.join(variable_dict['outdir'],variable_dict['run_name']+'_qc.csv')))
            printc("\n Pipeline complete\n", "HEADER")
            sys.exit(0)
        else:
            my_log.error("Something went wrong with the analysis. Investigate.")
            sys.exit(1)

##########
# %% run analysis                                                     #

if __name__ == '__main__':
    main()
