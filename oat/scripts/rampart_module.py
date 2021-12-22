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
import sys
import json
from datetime import date
import subprocess
import re
import webbrowser
import shlex

# =============================================================#
# %%GLOBAL VARIABLES AND DEPENDENCIES                            #
# =============================================================#
version = "0.7.2"
TODAY = date.today().strftime("%Y-%m-%d")
hash_line = "#" * 46

thisdir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.split(thisdir)[0]

# ============================================================#
# Define functions                                            #
# ============================================================#

# %% rampart_json
def rampart_json(variable_dict):
    """
    Prepare json files for RAMPART, one for each protocol
    """
    global my_log, run_data, run_name, rampart_outdir
    rampart_outdir = variable_dict["rampart_outdir"]
    my_log = variable_dict["my_log"]
    run_data = variable_dict["run_data"]
    run_name = variable_dict["run_name"]
    if not os.path.exists(rampart_outdir):
        os.makedirs(rampart_outdir)
    my_log.info(
        "Preparing json files for RAMPART"
    )
    protocol = variable_dict["protocol"]
    p_string = os.path.join(rampart_outdir, run_name + "_" + protocol)
    if os.path.exists(p_string):
        shutil.rmtree(p_string)
        os.makedirs(p_string)
    else:
        os.makedirs(p_string)
    samples = []
    for i in range(0, len(run_data)):
        if not variable_dict['demultiplexed']:
            tmp_dict = {
                "name": run_data["id"].iloc[i],
                "description": "",
                "barcodes": [run_data["barcode"].iloc[i]],
            }
        else:
            barcode = run_data["barcode"].iloc[i]
            fixed_barcode = re.sub("BC", "barcode", barcode)
            tmp_dict = {
                "name": run_data["id"].iloc[i],
                "description": "",
                "barcodes": [fixed_barcode],
            }
        samples.append(tmp_dict)
        rampart_input = {}
        rampart_input["title"] = run_name + "_" + protocol
        rampart_input["samples"] = samples
    fname = os.path.join(p_string, "run_configuration.json")
    with open(fname, "w") as fp:
        json.dump(rampart_input, fp, indent=2)
    return "RAMPART .json file(s) created"


# %% rampart_cmd
def rampart_run(variable_dict):
    """
    Run RAMPART, one process for each protocol
    """
    my_log.info("Running RAMPART")
    rampart_exe = shutil.which("rampart")
    protocol_name = variable_dict["protocol"]
    json_location = os.path.join(rampart_outdir, run_name + "_" + protocol_name)
    ports = "6010 6011"
    port_address = "http://localhost:6010"
    port_message = "RAMPART: view {0} protocol sequencing at {1}".format(
        protocol_name, port_address
    )
    webbrowser.open(port_address, new=1)
    my_log.info(port_message)
    protocol_path = os.path.join(script_dir, "protocols", protocol_name, "rampart")
    cmd = '{0} --verbose --protocol {1} --ports {2} --basecalledPath {3} --annotationOptions require_two_barcodes="False" barcode_set="rapid"'.format(
        rampart_exe, protocol_path, ports, variable_dict["basecalledPath"]
    )
    rampart_proc = subprocess.Popen(
        shlex.split(cmd),
        cwd=json_location,
        shell=False,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    fname = os.path.join(rampart_outdir, TODAY + "_" + run_name + "_RAMPART_cmd.txt")
    with open(fname, "w") as f:
        f.write("\n".join(map(str, cmd)))
    my_log.info("Opening default web browser")
    my_log.info("Sequencing information should appear after ~10 seconds")
    my_log.info("When finished monitoring the run with RAMPART, press enter")
    run_time = input("Press enter when ready to finish RAMPART ")
    rampart_proc.terminate()
    my_log.info("RAMPART commands written to: " + fname)
    my_log.info("RAMPART module complete")
    return ()
