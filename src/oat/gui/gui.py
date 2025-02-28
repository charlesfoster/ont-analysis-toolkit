#!/usr/bin/env python3
import sys
import os
import subprocess
import psutil
import re

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QGridLayout, QLineEdit,
    QLabel, QPushButton, QHBoxLayout, QGroupBox, QSpinBox, QDoubleSpinBox,
    QFileDialog, QCheckBox, QComboBox, QScrollArea, QMessageBox
)
from PyQt5.QtCore import Qt

def get_default_max_memory():
    """Return available memory in MB (rounded)."""
    return int(round(psutil.virtual_memory().available / (1024 * 1024)))

def generate_reference_info():
    """
    Dynamically generate a dictionary mapping reference IDs (derived from fasta filenames)
    to descriptions by scanning all .fasta files in the ../cli/references folder.
    """
    thisdir = os.path.abspath(os.path.dirname(__file__))
    refdir = os.path.join(thisdir, "..", "cli", "references")
    ref_dict = {}
    try:
        for filename in os.listdir(refdir):
            if filename.endswith(".fasta"):
                ref_id = filename.replace(".fasta", "")
                with open(os.path.join(refdir, filename), "r") as f:
                    header = f.readline().strip()
                    # Remove any trailing info after a comma
                    description = re.sub(r",.*", "", " ".join(header.split(" ")[1:]))
                ref_dict[ref_id] = description
    except Exception as e:
        # Fallback in case of error
        ref_dict = {"MN908947.3": "Default: SARS-CoV-2 (Wuhan Hu-1)"}
    return ref_dict

# Dictionaries for option descriptions
BARCODE_KIT_DESCRIPTIONS = {
    "SQK-RBK004": "Rapid 12-barcode kit",
    "SQK-RBK110-96": "Rapid 96-barcode kit",
    "SQK-RBK114-24": "Chemistry 14 rapid 12-barcode kit",
    "SQK-RBK114-96": "Chemistry 14 rapid 96-barcode kit",
    "EXP-NBD104": "Native 12-barcode kit",
    "EXP-NBD114": "Native 12-barcode expansion kit 13-24",
    "EXP-NBD196": "Native 96-barcode kit",
    "SQK-RPB004": "Rapid kit for adaptive sampling",
    "SQK-LSK114": "Ligation kit",
}

MODULE_DESCRIPTIONS = {
    "All": "Monitor run with RAMPART then run analysis",
    "Rampart": "Monitor run with RAMPART only",
    "Analysis": "Run analysis only",
}

GUPPY_MODEL_DESCRIPTIONS = {
    "r941_min_high_g360": "Chemistry v9, HAC basecalling",
    "r941_min_sup_g507": "Chemistry v9, SUP basecalling",
    "r1041_e82_260bps_sup_g632": "Chemistry v10, SUP basecalling",
    "Other": "Choose your own adventure: specify below",
}

CLAIR3_MODEL_DESCRIPTIONS = {
    "/opt/models/r941_prom_hac_g360+g422": "Default within clair3 container",
    "/opt/models/r941_prom_sup_g5014": "Super accurate model for v9.4.1 chemistry",
    "Other": "Choose your own adventure: specify below",
}

class AnalysisToolGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ONT Analysis Toolkit GUI v3")
        # Dynamically generate the reference info from fasta files.
        self.reference_dict = generate_reference_info()
        self.setupUI()
        self.resize(800, 550)

    def setupUI(self):
        # Create a scroll area for all content
        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        container = QWidget()
        scroll_area.setWidget(container)
        self.setCentralWidget(scroll_area)
        main_layout = QVBoxLayout(container)

        # ----- BASIC PARAMETERS GROUP -----
        basic_group = QGroupBox("Basic Parameters")
        basic_layout = QGridLayout()
        basic_group.setLayout(basic_layout)
        main_layout.addWidget(basic_group)
        row = 0

        # Samples File (CSV)
        basic_layout.addWidget(QLabel("Samples File (CSV):"), row, 0)
        self.samples_file = QLineEdit()
        samples_file_btn = QPushButton("Browse")
        samples_file_btn.clicked.connect(self.browseSamplesFile)
        samples_file_layout = QHBoxLayout()
        samples_file_layout.addWidget(self.samples_file)
        samples_file_layout.addWidget(samples_file_btn)
        basic_layout.addLayout(samples_file_layout, row, 1)
        basic_layout.addWidget(QLabel("CSV file with sample metadata."), row, 2)
        row += 1

        # Barcode Kit
        basic_layout.addWidget(QLabel("Barcode Kit:"), row, 0)
        self.combo_barcode_kit = QComboBox()
        self.combo_barcode_kit.addItems(list(BARCODE_KIT_DESCRIPTIONS.keys()))
        self.combo_barcode_kit.setCurrentText("SQK-RBK004")
        self.combo_barcode_kit.currentTextChanged.connect(
            lambda option: self.update_desc(self.label_barcode_desc, BARCODE_KIT_DESCRIPTIONS, option)
        )
        basic_layout.addWidget(self.combo_barcode_kit, row, 1)
        self.label_barcode_desc = QLabel(BARCODE_KIT_DESCRIPTIONS["SQK-RBK004"])
        basic_layout.addWidget(self.label_barcode_desc, row, 2)
        row += 1

        # Consensus Frequency
        basic_layout.addWidget(QLabel("Consensus Frequency:"), row, 0)
        self.consensus_freq = QDoubleSpinBox()
        self.consensus_freq.setDecimals(2)
        self.consensus_freq.setRange(0.0, 1.0)
        self.consensus_freq.setSingleStep(0.05)
        self.consensus_freq.setValue(0.0)
        basic_layout.addWidget(self.consensus_freq, row, 1)
        basic_layout.addWidget(QLabel("SNP incorporation threshold."), row, 2)
        row += 1

        # Indel Frequency
        basic_layout.addWidget(QLabel("Indel Frequency:"), row, 0)
        self.indel_freq = QDoubleSpinBox()
        self.indel_freq.setDecimals(2)
        self.indel_freq.setRange(0.0, 1.0)
        self.indel_freq.setSingleStep(0.05)
        self.indel_freq.setValue(0.40)
        basic_layout.addWidget(self.indel_freq, row, 1)
        basic_layout.addWidget(QLabel("Indel incorporation threshold."), row, 2)
        row += 1

        # Demultiplexed (checkbox)
        basic_layout.addWidget(QLabel("Demultiplexed:"), row, 0)
        self.demultiplexed = QCheckBox("Reads already demultiplexed")
        self.demultiplexed.setChecked(True)
        basic_layout.addWidget(self.demultiplexed, row, 1)
        basic_layout.addWidget(QLabel("Handled via MinKNOW."), row, 2)
        row += 1

        # Force Overwrite
        basic_layout.addWidget(QLabel("Force:"), row, 0)
        self.force = QCheckBox("Force overwrite of completed files")
        basic_layout.addWidget(self.force, row, 1)
        basic_layout.addWidget(QLabel("Overwrite previous outputs."), row, 2)
        row += 1

        # Dry Run
        basic_layout.addWidget(QLabel("Dry Run:"), row, 0)
        self.dry_run = QCheckBox("Dry run only")
        basic_layout.addWidget(self.dry_run, row, 1)
        basic_layout.addWidget(QLabel("Simulate analysis."), row, 2)
        row += 1

        # Module
        basic_layout.addWidget(QLabel("Module:"), row, 0)
        self.combo_module = QComboBox()
        self.combo_module.addItems(list(MODULE_DESCRIPTIONS.keys()))
        self.combo_module.setCurrentText("All")
        self.combo_module.currentTextChanged.connect(
            lambda option: self.update_desc(self.label_module_desc, MODULE_DESCRIPTIONS, option)
        )
        basic_layout.addWidget(self.combo_module, row, 1)
        self.label_module_desc = QLabel(MODULE_DESCRIPTIONS["All"])
        basic_layout.addWidget(self.label_module_desc, row, 2)
        row += 1

        # Threads
        basic_layout.addWidget(QLabel("Threads:"), row, 0)
        self.threads = QSpinBox()
        self.threads.setMinimum(1)
        self.threads.setMaximum(256)
        self.threads.setValue(psutil.cpu_count(logical=True))
        basic_layout.addWidget(self.threads, row, 1)
        basic_layout.addWidget(QLabel("CPU threads to use."), row, 2)
        row += 1

        # Reference Genome (dynamically generated from fasta files)
        basic_layout.addWidget(QLabel("Reference Genome:"), row, 0)
        self.reference = QComboBox()
        # Populate using keys from the dynamically generated dictionary
        self.reference.addItems(list(self.reference_dict.keys()))
        # Set default to "MN908947.3" if available
        if "MN908947.3" in self.reference_dict:
            index = self.reference.findText("MN908947.3")
            if index != -1:
                self.reference.setCurrentIndex(index)
        self.reference.currentTextChanged.connect(
            lambda option: self.update_desc(self.label_reference_desc, self.reference_dict, option)
        )
        basic_layout.addWidget(self.reference, row, 1)
        default_ref_desc = self.reference_dict.get(self.reference.currentText(), "")
        self.label_reference_desc = QLabel(default_ref_desc)
        basic_layout.addWidget(self.label_reference_desc, row, 2)
        row += 1

        # Variant Caller
        basic_layout.addWidget(QLabel("Variant Caller:"), row, 0)
        self.variant_caller = QComboBox()
        self.variant_caller.addItems(["clair3", "medaka"])
        self.variant_caller.setCurrentText("clair3")
        basic_layout.addWidget(self.variant_caller, row, 1)
        basic_layout.addWidget(QLabel("Select variant caller."), row, 2)
        row += 1

        # Minimum Depth
        basic_layout.addWidget(QLabel("Minimum Depth:"), row, 0)
        self.min_depth = QSpinBox()
        self.min_depth.setMinimum(1)
        self.min_depth.setMaximum(10000)
        self.min_depth.setValue(15)
        basic_layout.addWidget(self.min_depth, row, 1)
        basic_layout.addWidget(QLabel("Minimum read depth."), row, 2)
        row += 1

        # ----- ADVANCED OPTIONS GROUP -----
        # Add a separate checkbox to toggle the visibility of advanced options.
        self.advanced_toggle = QCheckBox("Show Advanced Options")
        self.advanced_toggle.toggled.connect(self.toggleAdvancedOptions)
        main_layout.addWidget(self.advanced_toggle)
        
        self.advanced_group = QGroupBox("Advanced Options")
        adv_layout = QGridLayout()
        self.advanced_group.setLayout(adv_layout)
        main_layout.addWidget(self.advanced_group)
        self.advanced_group.setVisible(False)  # Initially hidden
        arow = 0

        # Output Directory
        adv_layout.addWidget(QLabel("Output Directory:"), arow, 0)
        self.outdir = QLineEdit()
        outdir_btn = QPushButton("Browse")
        outdir_btn.clicked.connect(lambda: self.browseDirectory(self.outdir))
        outdir_layout = QHBoxLayout()
        outdir_layout.addWidget(self.outdir)
        outdir_layout.addWidget(outdir_btn)
        adv_layout.addLayout(outdir_layout, arow, 1)
        adv_layout.addWidget(QLabel("Where results are saved."), arow, 2)
        arow += 1

        # Rampart Outdir
        adv_layout.addWidget(QLabel("Rampart Outdir:"), arow, 0)
        self.rampart_outdir = QLineEdit(os.path.join(os.getcwd(), "rampart_files"))
        ramp_btn = QPushButton("Browse")
        ramp_btn.clicked.connect(lambda: self.browseDirectory(self.rampart_outdir))
        ramp_layout = QHBoxLayout()
        ramp_layout.addWidget(self.rampart_outdir)
        ramp_layout.addWidget(ramp_btn)
        adv_layout.addLayout(ramp_layout, arow, 1)
        adv_layout.addWidget(QLabel("Directory for Rampart files."), arow, 2)
        arow += 1

        # Print DAG
        adv_layout.addWidget(QLabel("Print DAG:"), arow, 0)
        self.print_dag = QCheckBox("Print workflow DAG")
        adv_layout.addWidget(self.print_dag, arow, 1)
        adv_layout.addWidget(QLabel("Display the DAG and exit."), arow, 2)
        arow += 1

        # Create Envs Only
        adv_layout.addWidget(QLabel("Create Envs:"), arow, 0)
        self.create_envs_only = QCheckBox("Create conda environments only")
        adv_layout.addWidget(self.create_envs_only, arow, 1)
        adv_layout.addWidget(QLabel("Only create environments, don't run analysis."), arow, 2)
        arow += 1

        # SNV Minimum Frequency
        adv_layout.addWidget(QLabel("SNV Minimum Frequency:"), arow, 0)
        self.snv_min_freq = QDoubleSpinBox()
        self.snv_min_freq.setDecimals(2)
        self.snv_min_freq.setRange(0.0, 1.0)
        self.snv_min_freq.setSingleStep(0.05)
        self.snv_min_freq.setValue(0.2)
        adv_layout.addWidget(self.snv_min_freq, arow, 1)
        adv_layout.addWidget(QLabel("Minimum allele frequency for SNVs."), arow, 2)
        arow += 1

        # Guppy Model (with custom option)
        adv_layout.addWidget(QLabel("Guppy Model:"), arow, 0)
        self.combo_guppy_model = QComboBox()
        self.combo_guppy_model.addItems(list(GUPPY_MODEL_DESCRIPTIONS.keys()))
        self.combo_guppy_model.setCurrentText("r941_min_high_g360")
        self.combo_guppy_model.currentTextChanged.connect(
            lambda opt: self.update_desc(self.label_guppy_desc, GUPPY_MODEL_DESCRIPTIONS, opt)
        )
        adv_layout.addWidget(self.combo_guppy_model, arow, 1)
        self.label_guppy_desc = QLabel(GUPPY_MODEL_DESCRIPTIONS["r941_min_high_g360"])
        adv_layout.addWidget(self.label_guppy_desc, arow, 2)
        arow += 1

        adv_layout.addWidget(QLabel("Custom Guppy Model:"), arow, 0)
        self.guypy_custom = QLineEdit()
        self.guypy_custom.setEnabled(False)
        adv_layout.addWidget(self.guypy_custom, arow, 1)
        self.combo_guppy_model.currentTextChanged.connect(self.updateGuppyCustom)
        arow += 1

        # Clair3 Model (with custom option)
        adv_layout.addWidget(QLabel("Clair3 Model:"), arow, 0)
        self.combo_clair3_model = QComboBox()
        self.combo_clair3_model.addItems(list(CLAIR3_MODEL_DESCRIPTIONS.keys()))
        self.combo_clair3_model.setCurrentText("/opt/models/r941_prom_hac_g360+g422")
        self.combo_clair3_model.currentTextChanged.connect(
            lambda opt: self.update_desc(self.label_clair3_desc, CLAIR3_MODEL_DESCRIPTIONS, opt)
        )
        adv_layout.addWidget(self.combo_clair3_model, arow, 1)
        self.label_clair3_desc = QLabel(CLAIR3_MODEL_DESCRIPTIONS["/opt/models/r941_prom_hac_g360+g422"])
        adv_layout.addWidget(self.label_clair3_desc, arow, 2)
        arow += 1

        adv_layout.addWidget(QLabel("Custom Clair3 Model:"), arow, 0)
        self.clair3_custom = QLineEdit()
        self.clair3_custom.setEnabled(False)
        adv_layout.addWidget(self.clair3_custom, arow, 1)
        self.combo_clair3_model.currentTextChanged.connect(self.updateClair3Custom)
        arow += 1

        # Alternate Analysis
        adv_layout.addWidget(QLabel("Alternate Analysis:"), arow, 0)
        self.alternate_analysis = QCheckBox("Alternate Analysis")
        adv_layout.addWidget(self.alternate_analysis, arow, 1)
        adv_layout.addWidget(QLabel("Enable alternate analysis mode."), arow, 2)
        arow += 1

        # Alt Cov Max
        adv_layout.addWidget(QLabel("Alt Cov Max:"), arow, 0)
        self.alt_cov_max = QDoubleSpinBox()
        self.alt_cov_max.setDecimals(2)
        self.alt_cov_max.setRange(0.0, 1000.0)
        self.alt_cov_max.setSingleStep(1.0)
        self.alt_cov_max.setValue(80.0)
        adv_layout.addWidget(self.alt_cov_max, arow, 1)
        adv_layout.addWidget(QLabel("Alternate analysis maximum coverage."), arow, 2)
        arow += 1

        # Alt Cov Min
        adv_layout.addWidget(QLabel("Alt Cov Min:"), arow, 0)
        self.alt_cov_min = QDoubleSpinBox()
        self.alt_cov_min.setDecimals(2)
        self.alt_cov_min.setRange(0.0, 1000.0)
        self.alt_cov_min.setSingleStep(1.0)
        self.alt_cov_min.setValue(40.0)
        adv_layout.addWidget(self.alt_cov_min, arow, 1)
        adv_layout.addWidget(QLabel("Alternate analysis minimum coverage."), arow, 2)
        arow += 1

        # Delete Reads
        adv_layout.addWidget(QLabel("Delete Reads:"), arow, 0)
        self.delete_reads = QCheckBox("Delete demultiplexed reads after analysis")
        adv_layout.addWidget(self.delete_reads, arow, 1)
        adv_layout.addWidget(QLabel("Remove reads post-analysis."), arow, 2)
        arow += 1

        # Redo Analysis
        adv_layout.addWidget(QLabel("Redo Analysis:"), arow, 0)
        self.redo_analysis = QCheckBox("Redo analysis (fresh run)")
        adv_layout.addWidget(self.redo_analysis, arow, 1)
        adv_layout.addWidget(QLabel("Delete output directory for a fresh run."), arow, 2)
        arow += 1

        # Additional nanoq parameters
        adv_layout.addWidget(QLabel("Additional nanoq parameters:"), arow, 0)
        self.additional_nanoq = QLineEdit()
        adv_layout.addWidget(self.additional_nanoq, arow, 1)
        adv_layout.addWidget(QLabel("Extra parameters for nanoq."), arow, 2)
        arow += 1

        # Skip Clipping
        adv_layout.addWidget(QLabel("Skip Clipping:"), arow, 0)
        self.skip_clipping = QCheckBox("Skip clipping of amplicon primers")
        adv_layout.addWidget(self.skip_clipping, arow, 1)
        adv_layout.addWidget(QLabel("Do not clip amplicon primers."), arow, 2)
        arow += 1

        # No Barcodes
        adv_layout.addWidget(QLabel("No Barcodes:"), arow, 0)
        self.no_barcodes = QCheckBox("No barcodes used during library prep")
        adv_layout.addWidget(self.no_barcodes, arow, 1)
        adv_layout.addWidget(QLabel("Specify if no barcodes were used."), arow, 2)
        arow += 1

        # MinKNOW Data Directory
        adv_layout.addWidget(QLabel("MinKNOW Data Directory:"), arow, 0)
        self.minknow_data = QLineEdit("/var/lib/minknow/data")
        minknow_btn = QPushButton("Browse")
        minknow_btn.clicked.connect(lambda: self.browseDirectory(self.minknow_data))
        minknow_layout = QHBoxLayout()
        minknow_layout.addWidget(self.minknow_data)
        minknow_layout.addWidget(minknow_btn)
        adv_layout.addLayout(minknow_layout, arow, 1)
        adv_layout.addWidget(QLabel("Directory for MinKNOW data."), arow, 2)
        arow += 1

        # No Update
        adv_layout.addWidget(QLabel("No Update:"), arow, 0)
        self.no_update = QCheckBox("Disable container version updating")
        adv_layout.addWidget(self.no_update, arow, 1)
        adv_layout.addWidget(QLabel("Prevent updating container versions."), arow, 2)
        arow += 1

        # List Protocols
        adv_layout.addWidget(QLabel("List Protocols:"), arow, 0)
        self.list_protocols = QCheckBox("List available protocols and exit")
        adv_layout.addWidget(self.list_protocols, arow, 1)
        adv_layout.addWidget(QLabel("Display protocols and exit."), arow, 2)
        arow += 1

        # Max Memory
        adv_layout.addWidget(QLabel("Max Memory (MB):"), arow, 0)
        self.max_memory = QSpinBox()
        self.max_memory.setMinimum(1)
        self.max_memory.setMaximum(100000)
        self.max_memory.setValue(get_default_max_memory())
        adv_layout.addWidget(self.max_memory, arow, 1)
        adv_layout.addWidget(QLabel("Maximum memory in MB."), arow, 2)
        arow += 1

        # Basecaller
        adv_layout.addWidget(QLabel("Basecaller:"), arow, 0)
        self.basecaller = QComboBox()
        self.basecaller.addItems(["dorado", "guppy"])
        self.basecaller.setCurrentText("dorado")
        adv_layout.addWidget(self.basecaller, arow, 1)
        adv_layout.addWidget(QLabel("Select the basecaller."), arow, 2)
        arow += 1

        # Rebasecall
        adv_layout.addWidget(QLabel("Rebasecall:"), arow, 0)
        self.rebasecall = QCheckBox("Rebasecall reads (for dorado)")
        adv_layout.addWidget(self.rebasecall, arow, 1)
        adv_layout.addWidget(QLabel("Enable rebasecalling."), arow, 2)
        arow += 1

        # Minimum QScore
        adv_layout.addWidget(QLabel("Minimum QScore:"), arow, 0)
        self.min_qscore = QSpinBox()
        self.min_qscore.setMinimum(1)
        self.min_qscore.setMaximum(100)
        self.min_qscore.setValue(9)
        adv_layout.addWidget(self.min_qscore, arow, 1)
        adv_layout.addWidget(QLabel("Minimum basecalling quality score."), arow, 2)
        arow += 1

        # Dorado Model
        adv_layout.addWidget(QLabel("Dorado Model (full path):"), arow, 0)
        self.dorado_model = QLineEdit()
        adv_layout.addWidget(self.dorado_model, arow, 1)
        adv_layout.addWidget(QLabel("Full path to the Dorado model."), arow, 2)
        arow += 1

        # Quiet
        adv_layout.addWidget(QLabel("Quiet:"), arow, 0)
        self.quiet = QCheckBox("Quiet mode (suppress commands)")
        adv_layout.addWidget(self.quiet, arow, 1)
        adv_layout.addWidget(QLabel("Suppress printing of commands."), arow, 2)
        arow += 1

        # Report
        adv_layout.addWidget(QLabel("Report:"), arow, 0)
        self.report = QCheckBox("Generate Snakemake report")
        adv_layout.addWidget(self.report, arow, 1)
        adv_layout.addWidget(QLabel("Generate a report after analysis."), arow, 2)
        arow += 1

        # ----- BUTTONS & STATUS -----
        btn_layout = QHBoxLayout()
        run_btn = QPushButton("Run Analysis")
        run_btn.setStyleSheet("background-color: green; color: white; font-weight: bold;")
        run_btn.clicked.connect(self.runAnalysis)
        btn_layout.addWidget(run_btn)
        reset_btn = QPushButton("Reset")
        reset_btn.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        reset_btn.clicked.connect(self.resetFields)
        btn_layout.addWidget(reset_btn)
        main_layout.addLayout(btn_layout)
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        main_layout.addWidget(self.status_label)

    # -----------------------
    # Helper Methods
    # -----------------------
    def update_desc(self, label, desc_dict, option):
        label.setText(desc_dict.get(option, ""))

    def browseSamplesFile(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Samples CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if file_path:
            self.samples_file.setText(file_path)

    def browseDirectory(self, target_line_edit):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            target_line_edit.setText(directory)

    def updateGuppyCustom(self, text):
        if text == "Other":
            self.guypy_custom.setEnabled(True)
        else:
            self.guypy_custom.setEnabled(False)
            self.guypy_custom.clear()

    def updateClair3Custom(self, text):
        if text == "Other":
            self.clair3_custom.setEnabled(True)
        else:
            self.clair3_custom.setEnabled(False)
            self.clair3_custom.clear()

    def toggleAdvancedOptions(self, checked):
        self.advanced_group.setVisible(checked)

    def generate_cli_input(self):
        """Generate a list of CLI arguments from the current GUI options."""
        args = []
        # Positional argument: samples file
        args.append(self.samples_file.text().strip())
        # Basic options
        args.extend(["--barcode_kit", self.combo_barcode_kit.currentText()])
        args.extend(["--consensus_freq", str(self.consensus_freq.value())])
        args.extend(["--indel_freq", str(self.indel_freq.value())])
        if self.demultiplexed.isChecked():
            args.append("--demultiplexed")
        if self.force.isChecked():
            args.append("--force")
        if self.dry_run.isChecked():
            args.append("--dry_run")
        args.extend(["--module", self.combo_module.currentText().upper()])
        args.extend(["--threads", str(self.threads.value())])
        args.extend(["--reference", self.reference.currentText().strip()])
        args.extend(["--variant_caller", self.variant_caller.currentText().strip()])
        args.extend(["--min_depth", str(self.min_depth.value())])
        # Advanced options (always added)
        if self.outdir.text().strip():
            args.extend(["--outdir", self.outdir.text().strip()])
        args.extend(["--rampart_outdir", self.rampart_outdir.text().strip()])
        if self.print_dag.isChecked():
            args.append("--print_dag")
        if self.create_envs_only.isChecked():
            args.append("--create_envs_only")
        args.extend(["--snv_min_freq", str(self.snv_min_freq.value())])
        guppy_model = (self.guypy_custom.text().strip() if self.combo_guppy_model.currentText() == "Other"
                       else self.combo_guppy_model.currentText())
        args.extend(["--guppy_model", guppy_model])
        clair3_model = (self.clair3_custom.text().strip() if self.combo_clair3_model.currentText() == "Other"
                        else self.combo_clair3_model.currentText())
        args.extend(["--clair3_model", clair3_model])
        if self.alternate_analysis.isChecked():
            args.append("--alternate_analysis")
        args.extend(["--alt_cov_max", str(self.alt_cov_max.value())])
        args.extend(["--alt_cov_min", str(self.alt_cov_min.value())])
        if self.delete_reads.isChecked():
            args.append("--delete_reads")
        if self.redo_analysis.isChecked():
            args.append("--redo_analysis")
        if self.additional_nanoq.text().strip():
            args.extend(["--additional_nanoq", self.additional_nanoq.text().strip()])
        if self.skip_clipping.isChecked():
            args.append("--skip_clipping")
        if self.no_barcodes.isChecked():
            args.append("--no_barcodes")
        if self.minknow_data.text().strip():
            args.extend(["--minknow_data", self.minknow_data.text().strip()])
        if self.no_update.isChecked():
            args.append("--no_update")
        if self.list_protocols.isChecked():
            args.append("--list_protocols")
        args.extend(["--max_memory", str(self.max_memory.value())])
        args.extend(["--basecaller", self.basecaller.currentText().strip()])
        if self.rebasecall.isChecked():
            args.append("--rebasecall")
        args.extend(["--min_qscore", str(self.min_qscore.value())])
        if self.dorado_model.text().strip():
            args.extend(["--dorado_model", self.dorado_model.text().strip()])
        if self.quiet.isChecked():
            args.append("--quiet")
        if self.report.isChecked():
            args.append("--report")
        return args

    def runAnalysis(self):
        # Validate mandatory field
        if not self.samples_file.text().strip():
            QMessageBox.critical(self, "Error", "Samples file is required.")
            return

        # Generate CLI input from GUI parameters
        cli_args = self.generate_cli_input()
        print("Command run:")
        print("oat" + " ".join(cli_args))

        # Close the GUI and invoke the CLI main() with the generated arguments
        self.close()
        import oat.cli.cli as cli
        sys.argv = [sys.argv[0]] + cli_args
        cli.main(sys.argv[1:])

    def resetFields(self):
        # Basic Options
        self.samples_file.clear()
        self.combo_barcode_kit.setCurrentText("SQK-RBK004")
        self.consensus_freq.setValue(0.0)
        self.indel_freq.setValue(0.40)
        self.demultiplexed.setChecked(True)
        self.force.setChecked(False)
        self.dry_run.setChecked(False)
        self.combo_module.setCurrentText("All")
        self.threads.setValue(psutil.cpu_count(logical=True))
        if "MN908947.3" in self.reference_dict:
            index = self.reference.findText("MN908947.3")
            if index != -1:
                self.reference.setCurrentIndex(index)
        self.variant_caller.setCurrentText("clair3")
        self.min_depth.setValue(15)
        # Advanced Options
        self.outdir.setText("")
        self.rampart_outdir.setText(os.path.join(os.getcwd(), "rampart_files"))
        self.print_dag.setChecked(False)
        self.create_envs_only.setChecked(False)
        self.snv_min_freq.setValue(0.2)
        self.combo_guppy_model.setCurrentText("r941_min_high_g360")
        self.guypy_custom.clear()
        self.guypy_custom.setEnabled(False)
        self.combo_clair3_model.setCurrentText("/opt/models/r941_prom_hac_g360+g422")
        self.clair3_custom.clear()
        self.clair3_custom.setEnabled(False)
        self.alternate_analysis.setChecked(False)
        self.alt_cov_max.setValue(80.0)
        self.alt_cov_min.setValue(40.0)
        self.delete_reads.setChecked(False)
        self.redo_analysis.setChecked(False)
        self.additional_nanoq.clear()
        self.skip_clipping.setChecked(False)
        self.no_barcodes.setChecked(False)
        self.minknow_data.setText("/var/lib/minknow/data")
        self.no_update.setChecked(False)
        self.list_protocols.setChecked(False)
        self.max_memory.setValue(get_default_max_memory())
        self.basecaller.setCurrentText("dorado")
        self.rebasecall.setChecked(False)
        self.min_qscore.setValue(9)
        self.dorado_model.clear()
        self.quiet.setChecked(False)
        self.report.setChecked(False)
        self.advanced_toggle.setChecked(False)  # This hides the advanced group
        self.status_label.setText("")

def main():
    app = QApplication(sys.argv)
    gui = AnalysisToolGUI()
    gui.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()