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

def check_directory_permissions(directory_path):
    if not os.access(directory_path, os.W_OK):
        print(f"The directory is write protected: {directory_path}")
        password = check_password()
        output = (True,password)
    else:
        output = (False,None)
    return output

def check_password():
    correct = False
    while not correct:
        password = getpass("Please enter your password: ")
        cmd = f" echo {password} | sudo -S echo ''"
        out, error = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).communicate()
        expected_error = f"[sudo] password for {os.environ['USER']}: "
        if error == expected_error:
            correct = True
        else:
            print("Incorrect password entered. Try again.")
    return password

def demultiplex_reads(variable_dict):
    global my_log, sample_dict, outdir, bcodeDir, protocol
    my_log = variable_dict["my_log"]
    run_data = variable_dict["run_data"]
    protocol = variable_dict["protocol"]

    my_log.info("Working with {0} protocol samples".format(protocol))
    protocol = protocol.upper()
    outdir = variable_dict["outdir"]
    sample_dict = dict(tuple(run_data.groupby("id")))
    bcodeDir = os.path.join(outdir, "barcodes")
    if not os.path.exists(bcodeDir):
        os.makedirs(bcodeDir)

    basecall_tool = variable_dict['basecaller']
    guppy_barcoder = shutil.which("guppy_barcoder")
    # dorado options
    dorado = shutil.which("dorado")
    if dorado is None:
        if guppy_barcoder is None:
            guppy_barcoder = shutil.which("ont_barcoder")
            if guppy_barcoder is None:
                my_log.error("No barcoder detected: exiting")
                sys.exit()
    basecalling_needed = variable_dict["rebasecall"]
    min_qscore = variable_dict["min_qscore"]
    basecall_model = variable_dict["dorado_model"]
    basecalledPath = variable_dict["basecalledPath"]
    run_dir = os.path.split(basecalledPath)[0]
    pod5_dir = os.path.join(os.path.split(basecalledPath)[0],"pod5_pass")
    if not os.path.exists(pod5_dir):
        os.makedirs(pod5_dir)
    oat_analyses = os.path.join(run_dir,'oat_analyses')
    if not os.path.exists(oat_analyses):
        os.makedirs(oat_analyses)
    # oat_basecall_dir = os.path.join(oat_analyses, "oat_basecalling")
    oat_basecall_dir = os.path.join(outdir, "oat_basecalling")
    if not os.path.exists(oat_basecall_dir):
        os.makedirs(oat_basecall_dir)
    basecalled_bam = os.path.join(oat_basecall_dir,"oat_basecalled.bam")
    # demux_outdir = os.path.join(oat_analyses, "oat_demux")
    demux_outdir = bcodeDir
    if not os.path.exists(demux_outdir):
        os.makedirs(demux_outdir)
    barcode_kit_name = variable_dict["barcode_kit"]
    if "SQK-RBK" in barcode_kit_name:
        barcode_type = "rapid"
        barcode_option = ""
    elif "EXP-NBD" in barcode_kit_name:
        barcode_type = "ligation"
        barcode_option = "--require_barcodes_both_ends"
    else:
        my_log.warning("Unsure of barcode_kit type - assuming barcodes not required at both ends")
        barcode_type = "rapid"
        barcode_option = ""
    # check for fastqs in bcodeDir already
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
            return(variable_dict)
        elif decision == 2:
            # remove the need for sudo everything
            password_check = check_directory_permissions(run_dir)
            if password_check[0]:
                print(f"Temporarily granting {os.environ['USER']} permission to write in {run_dir}...")
                cmd = f"echo {password_check[1]} | sudo --stdin chown {os.environ['USER']}:{os.environ['USER']} {run_dir}"
                os.system(cmd)
            # my_log.info("Temporarily stopping MinKNOW service to free up GPU memory")
            try:
                password
            except:
                # password = getpass("Please enter your password: ")
                password = password_check[1]
                variable_dict["password"] = password
                # cmd = "echo {} | sudo --stdin service minknow stop".format(password)
                # subprocess.Popen(
                #     shlex.split(cmd),
                #     shell=False,
                #     stderr=subprocess.PIPE,
                #     stdout=subprocess.PIPE,
                #     encoding="ascii",
                # )
                # my_log.info("Sleeping for 10 seconds to let the changes kick in")
                # time.sleep(10)
            try:
                if basecall_tool == 'dorado':
                    if basecalling_needed:
                        my_log.info("Basecalling with dorado.")
                        cmd = "{0} basecaller {1} {2} --recursive --kit-name {3} --min-qscore {4} {5} > {6}".format(
                                dorado, basecall_model, pod5_dir, barcode_kit_name, min_qscore, barcode_option, basecalled_bam
                            )
                        my_log.debug(cmd)
                        subprocess.Popen(
                            # shlex.split(cmd),
                            # shell=False,
                            cmd,
                            shell=True,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            universal_newlines=True,
                        ).communicate()
                        demux_cmd = "{0} demux --output-dir {1} --emit-fastq --no-classify {2}".format(
                            dorado, demux_outdir, basecalled_bam
                        )
                    else:
                        # make sure there is a monolithic pod5 file
                        pod5s = [x for x in os.listdir(pod5_dir) if x.endswith(".pod5")]
                        if len(pod5s) > 1:
                            my_log.info("Merging pod5 files.")
                            merged_pod5 = os.path.join(demux_outdir,"oat_merged.pod5")
                            cmd = "pod5 merge {0} -o {1}".format(
                                os.path.join(pod5_dir,"*.pod5"), merged_pod5
                            )
                        else:
                            merged_pod5 = pod5s[0]
                        demux_cmd = "{0} demux --output-dir {1} --emit-fastq --no-classify {2}".format(
                            dorado, demux_outdir, merged_pod5
                        )
                    my_log.info("Demultiplexing with dorado.")
                    my_log.debug(demux_cmd)
                    subprocess.Popen(
                        # shlex.split(demux_cmd),
                        # shell=False,
                        demux_cmd,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        universal_newlines=True,
                    ).communicate()
                    # reformat outfiles into their own directories
                    fastq_files = [x for x in os.listdir(demux_outdir) if x.endswith(".fastq")]
                    for file in fastq_files:
                        barcode = file.split('_')[-1].replace('.fastq', '')
                        dir_path = os.path.join(demux_outdir, barcode)
                        if not os.path.exists(dir_path):
                            os.makedirs(dir_path)
                        original_file_path = os.path.join(demux_outdir, file)
                        new_file_path = os.path.join(dir_path, file)
                        shutil.move(original_file_path, new_file_path)                        
                    variable_dict['basecalledPath'] = demux_outdir
                    print(variable_dict['basecalledPath'])
                elif basecall_tool == 'guppy':
                    my_log.info("Demultiplexing with guppy_barcoder.")
                    cmd = "{0} --input_path {1} --save_path {2} --detect_mid_strand_barcodes -x auto --barcode_kits {3}  {4}".format(
                        guppy_barcoder, basecalledPath, bcodeDir, barcode_kit_name, barcode_option
                    )
                else:
                    my_log.error("Invalid basecall tool: only 'dorado' or 'guppy' supported")
                    sys.exit(1)
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
                # # just in case gpu can't be detected
                # cmd = "{0} --input_path {1} --save_path {2} --detect_mid_strand_barcodes --barcode_kits {3} --trim_barcodes {4}".format(
                #     guppy_barcoder, basecalledPath, bcodeDir, barcode_kit_name, barcode_option
                # )
                # my_log.debug(cmd)
                # subprocess.Popen(
                #     shlex.split(cmd),
                #     shell=False,
                #     stdin=subprocess.PIPE,
                #     stdout=subprocess.PIPE,
                #     stderr=subprocess.STDOUT,
                #     universal_newlines=True,
                # ).communicate()
                my_log.error("Something went wrong. Is your GPU working? No CUDA errors?")
                sys.exit(1)
            barcode_dirs = next(os.walk(bcodeDir))[1]
            unwanted_dirs = [
                x for x in barcode_dirs if x not in variable_dict["barcodes_used"]
            ]
            for bdir in unwanted_dirs:
                shutil.rmtree(os.path.join(bcodeDir, bdir))
            # os.remove(os.path.join(bcodeDir, "barcoding_summary.txt"))
            # restoring ROOT permission if necessary
            if password_check[0]:
                print(f"\nChanging directory permissions back to ROOT...")
                cmd = f"echo {password_check[1]} | sudo --stdin chown root:root {run_dir}"
                os.system(cmd)
            my_log.info("Demultiplexing complete")
            return(variable_dict)
        elif decision == 3:
            my_log.error("Quitting")
            sys.exit()
        else:
            my_log.error("The only valid options are 1 or 2. Try again.")


def filter_reads(variable_dict):
    my_log.info("Filtering reads by length with 'nanoq'.")
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
                cmd = "nanoq --min-len {0} --max-len {1} -i {2} >> {3}".format(
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
    my_log = variable_dict["my_log"]
    outdir = variable_dict["outdir"]
    variable_dict["reads_dir"] = os.path.join(outdir, "reads")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(os.path.join(outdir, "reads")):
        os.makedirs(os.path.join(outdir, "reads"))
    run_data = variable_dict["run_data"]
    bcodeDir = os.path.join(outdir, "barcodes")
    if not os.path.exists(bcodeDir):
        os.makedirs(bcodeDir)
    sample_dict = dict(tuple(run_data.groupby("id")))
    if variable_dict['demultiplexed']:
        my_log.info("Option 'demultiplexed' chosen")
        my_log.info("Checking to see if reads previously demultiplexed and filtered by 'oat'")
        outfile_list = []
        for sample in sample_dict:
            outfile_list.append(os.path.join(outdir, "reads", sample + ".fastq"))
        if all([os.path.exists(x) for x in outfile_list]):
            my_log.info("All reads found. Continuing.")
            return
        else:
            my_log.warning("Not all reads can be located in the 'oat' outdir")
            my_log.warning("Working with reads in {0}".format(variable_dict['basecalledPath']))
    variable_dict["my_log"].info("Filtering reads by length with 'nanoq'")
    for sample in sample_dict:
        bcode = "".join(sample_dict[sample]["barcode"].to_list()).replace(
            "BC", "barcode"
        )
        fqdir = os.path.join(bcodeDir, bcode)
        if not os.path.isdir(fqdir):
            os.makedirs(fqdir)
        outfile = os.path.join(outdir, "reads", sample + ".fastq")
        if not os.path.isfile(outfile):
            readsDir = os.path.join(variable_dict["basecalledPath"] + "/" + bcode)
            print(readsDir)
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
            for fastq in fastqs:
                print(fastq)
                cmd = "nanoq --min-len {0} --max-len {1} -i {2} >> {3}".format(
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
        else:
            f"{outfile} already exists - skipping filtering..."
