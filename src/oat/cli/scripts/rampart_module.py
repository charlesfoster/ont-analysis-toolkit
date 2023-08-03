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
    if variable_dict['no_barcodes']:
        extra_options = 'barcode_threshold=100 require_two_barcodes="True" limit_barcodes_to="BC01' # ensuring everything becomes unassigned
    else:
        extra_options = 'require_two_barcodes="False"'
    cmd = '{0} --verbose --protocol {1} --ports {2} --basecalledPath {3} --annotationOptions barcode_set="rapid" {4}'.format(
        rampart_exe, protocol_path, ports, variable_dict["basecalledPath"], extra_options
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
       # f.write("\n".join(map(str, cmd)))
        f.write(cmd)
    my_log.info("Opening default web browser")
    my_log.info("Sequencing information should appear after ~10 seconds")
    my_log.info("When finished monitoring the run with RAMPART, press enter")
    run_time = input("Press enter when ready to finish RAMPART ")
    rampart_proc.terminate()
    my_log.info("RAMPART commands written to: " + fname)
    my_log.info("RAMPART module complete")
    return ()

# %% To handle the case when adaptive sampling is run without barcoding
# modified from advice here: https://medium.com/@sairajtri/magic-of-moving-the-files-among-folders-with-python-using-watchdog-f0dadfbb7a07

import time
import os
import sys
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from getpass import getpass

class movehandler(FileSystemEventHandler):
    # take password tuple as argument
    def __init__(self, password_check, pathtowatch, destination):
        self.password_check = password_check
        self.pathtowatch = pathtowatch
        self.destination = destination
    
    #overriding the on_modified method
    def on_modified(self, event):
        dir1 = [x for x in os.listdir(self.pathtowatch) if x.endswith(".fastq")]
        for fname in dir1:
            fname1=fname
            fname2=fname
            source = os.path.join(self.pathtowatch,fname1)
            newdestination = os.path.join(self.destination,fname2)
            file_done = False
            file_size = -1
            while file_size != os.path.getsize(self.pathtowatch):
                file_size = os.path.getsize(self.destination)
                time.sleep(1)
            while not file_done:
                try:
                    print(f"Read found: {fname} ")
                    print(f"Moving to barcode folder: {self.destination} ")
                    if self.password_check[0]:
                        cmd = f"echo {self.password_check[1]} | sudo --stdin mv {source} {newdestination}"
                    else:
                        cmd = f"mv {source} {newdestination}"
                    os.system(cmd)
                    file_done = True
                except:
                    return True        

def check_directory_permissions(directory_path):
    if not os.access(directory_path, os.W_OK):
        print(f"The directory is write protected: {directory_path}")
        password = getpass("Please enter your password: ")
        output = (True,password)
    else:
        output = (False,None)
    return output

def initiate_watchdog(pathtowatch,destination,password_check):
    #creating event handler object
    event_handler = movehandler(password_check,pathtowatch,destination)
    #creating observer object
    observer = Observer()
    observer.schedule(event_handler, pathtowatch, recursive=True)
    observer.start()
    try:
        while True:
            print("Watching for reads: press CTRL+C to end the watch and end RAMPART")
            time.sleep(10)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

# %%
def rampart_watchdog(variable_dict):
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
    rampart_cmd = '{0} --verbose --protocol {1} --ports {2} --basecalledPath {3} --annotationOptions barcode_set="rapid" require_two_barcodes="False"'.format(
        rampart_exe, protocol_path, ports, variable_dict["basecalledPath"]
    )
    fname = os.path.join(rampart_outdir, TODAY + "_" + run_name + "_RAMPART_cmd.txt")
    with open(fname, "w") as f:
        f.write(rampart_cmd)
    my_log.info("Opening default web browser")
    my_log.info("Sequencing information should appear after ~10 seconds")
    if not variable_dict['no_barcodes']:
        # start rampart proc
        rampart_proc = subprocess.Popen(
            shlex.split(cmd),
            cwd=json_location,
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        my_log.info("When finished monitoring the run with RAMPART, press enter")
        run_time = input("Press enter when ready to finish RAMPART ")
        rampart_proc.terminate()
    else:
        pathtowatch = variable_dict["basecalledPath"]
        destination = os.path.join(pathtowatch, variable_dict["pseudo_barcode"]) # need to create this elsewhere
        password_check = check_directory_permissions(pathtowatch)
        if not os.path.exists(destination):
            if password_check[0]:
                cmd = f"echo {password_check[1]} | sudo mkdir -p {destination}"
                os.system(cmd)
            else:
                os.makedirs(destination)
        # start rampart proc
        rampart_proc = subprocess.Popen(
            shlex.split(rampart_cmd),
            cwd=json_location,
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # start the watchdog process
        initiate_watchdog(pathtowatch,destination,password_check)
        rampart_proc.terminate()
    my_log.info("RAMPART commands written to: " + fname)
    my_log.info("RAMPART module complete")
    return ()
