#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last edited on Thur December 02 2021

@author: Dr Charles Foster
"""

# ==============================================================#
# %% Import Modules                                               #
# ==============================================================#
import os
import shutil
from datetime import date
import subprocess
import pathlib
import shlex
from getpass import getpass
import time
import glob
import sys

# =============================================================#
# %%GLOBAL VARIABLES AND DEPENDENCIES                            #
# =============================================================#
version = "0.7.2"
TODAY = date.today().strftime("%Y-%m-%d")
hash_line = "#" * 46

thisdir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.split(thisdir)[0]
main_dir = pathlib.Path(__file__).parent.parent.parent.resolve()


# ============================================================#
# Define functions                                            #
# ============================================================#

# %% artic_analysis - calls several of the above functions
def demultiplex_reads(variable_dict):
    global my_log, sample_dict, outdir, bcodeDir, protocol
    my_log = variable_dict["my_log"]
    run_data = variable_dict["run_data"]
    protocol = variable_dict["protocol"]
    guppy_barcoder = shutil.which("guppy_barcoder")
    basecalledPath = variable_dict["basecalledPath"]
    barcode_kit_name = variable_dict["barcode_kit_name"]

    my_log.info("Working with {0} protocol samples".format(protocol))
    protocol = protocol.upper()
    outdir = variable_dict["outdir"]
    sample_dict = dict(tuple(run_data.groupby("id")))
    bcodeDir = os.path.join(outdir, "barcodes")
    if not os.path.exists(bcodeDir):
        os.makedirs(bcodeDir)
    # check for fastqs in there already
    fq_check = glob.glob(bcodeDir + "/**/*fastq")
    decision = int()
    while decision not in [1, 2, 3]:
        if len(fq_check) > 0:
            my_log.warning(
                "Demultiplexing appears to have (at least partially) occurred"
            )
            try:
                decision = int(
                    input(
                        "\nDo you want to:\n(1) assume demultiplexing finished successfully\n(2) run demultiplexing from scratch\n(3) quit\n\nType 1, 2 or 3 then press enter: "
                    )
                )
            except:
                decision = 0
        else:
            decision = 2
        if decision == 1:
            my_log.info("Reads already demultiplexed. Skipping.")
            if os.path.isfile(os.path.join(bcodeDir, "barcoding_summary.txt")):
                os.remove(os.path.join(bcodeDir, "barcoding_summary.txt"))
        elif decision == 2:
            my_log.info("Temporarily stopping MinKNOW service to free up GPU memory")
            try:
                password
            except:
                password = getpass("Please enter your password: ")
                variable_dict["password"] = password
                cmd = "echo {} | sudo --stdin service minknow stop".format(password)
                subprocess.Popen(
                    shlex.split(cmd),
                    shell=False,
                    stderr=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    encoding="ascii",
                )
                my_log.info("Sleeping for 10 seconds to let the changes kick in")
                time.sleep(10)
            my_log.info("Demultiplexing with guppy_barcoder.")
            try:
                cmd = "{0} --input_path {1} --save_path {2} --detect_mid_strand_barcodes -x auto --barcode_kits {3}".format(
                    guppy_barcoder, basecalledPath, bcodeDir, barcode_kit_name
                )
                my_log.debug(cmd)
                subprocess.Popen(
                    shlex.split(cmd),
                    shell=False,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                ).communicate()
            except:
                # just in case gpu can't be detected
                cmd = "{0} --input_path {1} --save_path {2} --detect_mid_strand_barcodes --barcode_kits {3}".format(
                    guppy_barcoder, basecalledPath, bcodeDir, barcode_kit_name
                )
                my_log.debug(cmd)
                subprocess.Popen(
                    shlex.split(cmd),
                    shell=False,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                ).communicate()
            barcode_dirs = next(os.walk(bcodeDir))[1]
            unwanted_dirs = [
                x for x in barcode_dirs if x not in variable_dict["barcodes_used"]
            ]
            for bdir in unwanted_dirs:
                shutil.rmtree(os.path.join(bcodeDir, bdir))
            os.remove(os.path.join(bcodeDir, "barcoding_summary.txt"))
            my_log.info("Demultiplexing complete")
        elif decision == 3:
            my_log.error("Quitting")
            sys.exit()
        else:
            my_log.error("The only valid options are 1 or 2. Try again.")


def filter_reads(variable_dict):
    my_log.info("Filtering reads by length with 'filtlong'.")
    if not os.path.exists(os.path.join(outdir, "reads")):
        os.makedirs(os.path.join(outdir, "reads"))
    for sample in sample_dict:
        bcode = "".join(sample_dict[sample]["barcode"].to_list()).replace(
            "BC", "barcode"
        )
        fqdir = os.path.join(bcodeDir, bcode)
        outfile = os.path.join(outdir, "reads", sample + ".fastq")
        fastqs = glob.glob(fqdir + "/*.fastq")
        if os.path.exists(outfile):
            my_log.warning(
                f"Filtered reads file for {sample} already exists ({outfile})"
            )
            my_log.warning(
                "Please remove the file and try again if you want to generate the file fresh"
            )
        else:
            for fastq in fastqs:
                cmd = "filtlong --min_length {0} --max_length {1} {2} >> {3}".format(
                    variable_dict["min_len"], variable_dict["max_len"], fastq, outfile
                )
                subprocess.Popen(
                    cmd,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                ).communicate()
                # create symlinks for snakemake analysis
                cmd = "ln -s {0} {1}".format(
                    outfile, os.path.join(bcodeDir, bcode, os.path.basename(outfile))
                )
                subprocess.Popen(
                    shlex.split(cmd),
                    shell=False,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                ).communicate()
            # clean up unfiltered reads
        [os.remove(os.path.join(fqdir, x)) for x in os.listdir(fqdir) if "_runid_" in x]
    variable_dict["reads_dir"] = os.path.join(outdir, "reads")


def relocate_and_filter_reads(variable_dict):
    outdir = variable_dict["outdir"]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    run_data = variable_dict["run_data"]
    bcodeDir = os.path.join(outdir, "barcodes")
    if not os.path.exists(bcodeDir):
        os.makedirs(bcodeDir)
    sample_dict = dict(tuple(run_data.groupby("id")))
    for sample in sample_dict:
        bcode = "".join(sample_dict[sample]["barcode"].to_list()).replace(
            "BC", "barcode"
        )
        fqdir = os.path.join(bcodeDir, bcode)
        if not os.path.isdir(fqdir):
            os.makedirs(fqdir)
        if not os.path.isfile(os.path.join(outdir, sample + ".fastq")):
            readsDir = os.path.join(variable_dict["basecalledPath"] + "/" + bcode)
            if not os.path.exists(readsDir):
                variable_dict["my_log"].error(
                    "Cannot find demultiplexed reads for {0} in the expected location ({1})".format(
                        sample, readsDir
                    )
                )
                variable_dict["my_log"].error(
                    "Ensure that demultiplexing worked successfully with MinKNOW, or demultiplex again using this pipeline with option '-d'"
                )
                sys.exit(1)
            fastqs = [
                readsDir + "/" + file
                for file in os.listdir(readsDir)
                if file.endswith(".fastq")
            ]
            outfile = os.path.join(outdir, sample + ".fastq")
            for fastq in fastqs:
                cmd = "filtlong --min_length {0} --max_length {1} {2} >> {3}".format(
                    variable_dict["min_len"], variable_dict["max_len"], fastq, outfile
                )
                subprocess.Popen(
                    cmd,
                    shell=False,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                ).communicate()
            cmd = "ln -s {0} {1}".format(
                outfile, os.path.join(fqdir, os.path.basename(outfile))
            )
            subprocess.Popen(
                shlex.split(cmd),
                shell=False,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            ).communicate()
    variable_dict["reads_dir"] = os.path.join(outdir, "reads")
