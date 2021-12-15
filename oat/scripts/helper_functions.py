#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 10:31:05 2021

@author: cfos
"""
# import configparser
import sys
import os
import logging
from colorlog import ColoredFormatter
import pandas as pd
import re
import shutil
from datetime import date
import pathlib
import glob

TODAY = date.today().strftime("%Y-%m-%d")
thisdir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.split(thisdir)[0]
main_dir = pathlib.Path(__file__).parent.parent.parent.resolve()
reference_dir = os.path.join(
    pathlib.Path(__file__).parent.parent.resolve(), "references"
)

# %% printc
def printc(thing, level):
    """
    Print in colour :)
    """
    cols = {"HEADER": "\033[95m", "STD": "\033[96m"}
    col = cols[level]
    print(f"{col}{thing}\033[0m")
    return ()


# %% initiate_colorlog()
def initiate_colorlog(variable_dict, directory):
    global my_log, logfile
    """
    Set up logging to file with colored messages
    """
    logfile = directory + "/" + TODAY + "_analysis.log"
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M",
        filename=logfile,
        filemode="w",
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    LOGFORMAT = "  %(log_color)s%(asctime)s %(levelname)-8s%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(
        LOGFORMAT,
        log_colors={
            "DEBUG": "cyan",
            "INFO": "green",
            "WARNING": "yellow",
            "ERROR": "red",
            "CRITICAL": "red,bg_white",
        },
        datefmt="%Y-%m-%d %H:%M",
    )
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    my_log = logging.getLogger("ONT Pipeline")
    variable_dict["my_log"] = my_log
    variable_dict["logfile"] = logfile
    return ()


# %% check parameters
def check_arguments(variable_dict, args):
    """
    Check to make sure command line arguments are valid
    """
    global run_data, run_name, neg_controls, barcodes_used, barcode_kit_name

    if shutil.which("rampart") is None:
        my_log.error(
            "Some necessary programs cannot be detected. Have you activated the conda environment?"
        )
        sys.exit()
    if not os.path.isfile("".join(args.samples_file)):
        my_log.error(
            "Input file does not exist. Please check the filename and try again."
        )
        sys.exit()
    else:
        run_data = pd.read_csv("".join(args.samples_file), sep=",")
        run_data["id"] = run_data["id"].apply(str)
        if any(bool(re.search(r"\s", item)) for item in run_data["id"]):
            my_log.error(
                "Sample IDs cannot have spaces in their names. Please rename and try again. Offending sample ID(s):"
            )
            [
                print("\t" + item + "\n")
                for item in run_data["id"]
                if bool(re.search(r"\s", item))
            ]
            sys.exit()
        run_data["protocol"] = run_data["protocol"].str.upper()
        variable_dict["run_data"] = run_data
        if len(list(run_data["protocol"].apply(str).unique())) != 1:
            my_log.error(
                "A protocol must be specified in your input spreadsheet, and only one protocol is allowed."
            )
            sys.exit()
        protocol = "".join(list(run_data["protocol"].apply(str).unique()))
        variable_dict["protocol"] = protocol
        protocols_dir = os.path.join(
            pathlib.Path(__file__).parent.parent.resolve(), "protocols"
        )
        if not os.path.exists(os.path.join(protocols_dir, protocol)):
            my_log.error(
                "The specified protocol does not appear to exist. Please install the protocol appropriately within {}".format(
                    protocols_dir
                )
            )
            sys.exit()
        variable_dict["protocol_bed"] = glob.glob(
            os.path.join(protocols_dir, protocol) + "/**/*.scheme.bed", recursive=True
        )
        run_name = "".join(list(run_data["run_name"].apply(str).unique()))
        variable_dict["run_name"] = run_name
        neg_controls = run_data["id"][run_data["neg_control"] == True].tolist()
        variable_dict["neg_controls"] = neg_controls
        barcodes_used = [
            x.replace("BC", "barcode") for x in run_data["barcode"].to_list()
        ]
        variable_dict["barcodes_used"] = barcodes_used
    if not "".join(args.module) in ["RAMPART", "ANALYSIS", "ALL"]:
        my_log.error(
            "Module does not exist. Allowed values are 'rampart', 'analysis', or 'all'. Please check the spelling and try again."
        )
        sys.exit()
    if "".join(args.barcode_kit) == "12":
        barcode_kit_name = "SQK-RBK004"
    elif "".join(args.barcode_kit) == "96":
        barcode_kit_name = "SQK-RBK110-96"
    else:
        my_log.error(
            "Could not determine barcode kit that you used (only '12' or '96' supported as options)"
        )
        sys.exit()
    variable_dict["barcode_kit_name"] = barcode_kit_name
    if args.demultiplex:
        variable_dict["demultiplex"] = True
    else:
        variable_dict["demultiplex"] = False
    #    if not args.outdir:
    #        variable_dict['outdir'] = os.path.join(main_dir, "analysis_results", variable_dict['run_name'])

    reference = os.path.join(reference_dir, variable_dict["reference"] + ".fasta")
    annotation = os.path.join(reference_dir, variable_dict["reference"] + ".gff3")
    if not os.path.isfile(reference) or not os.path.isfile(annotation):
        my_log.error(
            "Reference files (fasta; gff3) do not appear to exist in the correct location. Please check and try again."
        )
        sys.exit()
    variable_dict["reference"] = reference
    variable_dict["annotation"] = annotation
    my_log.info("Run metadata: " + "".join(args.samples_file))
    my_log.info("Running the " + "".join(args.module + " module"))
    my_log.debug("Command used: " + " ".join(sys.argv))
    files = {"Run data": args.samples_file, "Module": "".join(args.module)}
    my_log.debug(files)

    return "Command line arguments are valid. Proceeding."


# %% find_runDirs
def find_runDirs(variable_dict, script_dir, minknow_dir):
    """
    Find directories associated with a sequencing run
    """
    global basecalledPath, fast5Dir
    # remove annotations directory, if it exists
    if os.path.isdir(os.path.join(script_dir, "annotations")):
        shutil.rmtree(os.path.join(script_dir, "annotations"))
    for dirpath, subdirs, files in os.walk(minknow_dir):
        for x in subdirs:
            if x.endswith("fastq_pass") and "/" + run_name + "/" in os.path.join(
                dirpath
            ):
                basecalledPath = os.path.join(dirpath, x)
                variable_dict["basecalledPath"] = basecalledPath
            if x.endswith("fast5_pass") and run_name + "/" in os.path.join(dirpath):
                fast5Dir = os.path.join(dirpath, x)
                variable_dict["fast5Dir"] = fast5Dir
    # if 'basecalledPath' not in locals() and 'fast5Dir' not in locals():
    if "basecalledPath" not in globals():
        msg = 'fastq_pass/fast5_pass not found for MinKNOW run {0}. Wait and try again, or check the "minknow_directory" in your config file.'.format(
            run_name
        )
        my_log.error(msg)
        sys.exit(1)
    return ()
