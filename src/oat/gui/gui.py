from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QGridLayout, QCheckBox,
                             QLabel, QComboBox, QLineEdit, QPushButton, QFileDialog, QMessageBox, QGroupBox)
import sys
import psutil
import os
from functools import partial

## FUNCTIONS
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

def check_env_variables():
    gmodel_check = os.getenv('GUPPY_MODEL')
    if gmodel_check is None:
        gmodel_check = 'r941_min_high_g360'
    c3model_check = os.getenv('CLAIR3_MODEL')
    if c3model_check is None:
        c3model_check = '/opt/models/r941_prom_hac_g360+g422'
    return([gmodel_check,c3model_check])

## DEFINE VARIABLES ##
important_models = check_env_variables()

BARCODE_KIT_DESCRIPTIONS = {
    "SQK-RBK004": "Rapid 12-barcode kit",
    "SQK-RBK110-96": "Rapid 96-barcode kit",
    "SQK-RBK114-24": "Chemistry 14 rapid 12-barcode kit",
    "SQK-RBK114-96": "Chemistry 14 rapid 96-barcode kit",
    "EXP-NBD104": "Native 12-barcode kit",
    "EXP-NBD114": "Native 12-barcode expansion kit 13-24",
    "EXP-NBD196": "Native 96-barcode kit",
}

DEMULTIPLEXED_DESCRIPTIONS = {
    "True": "Sample demultiplexing handled via MinKNOW",
    "False": "Sample demultiplexing not handled via MinKNOW",
}

SKIP_CLIPPING_DESCRIPTIONS = {
    "True": "Do not attempt to clip amplicon primers (useful for capture data)",
    "False": "Do not skip clipping of amplicon primers (recommended for amplicon data)",
}

MODULE_DESCRIPTIONS = {
    "All": "Monitor run with RAMPART then run analysis",
    "Rampart": "Monitor run with RAMPART only",
    "Analysis": "Run analysis only",
}

REFERENCE_DESCRIPTIONS = {
    "MN908947.3": "SARS-CoV-2 (Wuhan Hu-1)",
    "NC_006273.2": "Human cytomegalovirus (Merlin)",
}

GUPPY_MODEL_DESCRIPTIONS = {
    "r941_min_high_g360": "Chemistry version 9, HAC basecalling",
    "r941_min_sup_g507": "Chemistry version 9, SUP basecalling",
    "Other": "Choose your own adventure: specify below",
}

if os.getenv('GUPPY_MODEL') is not None:
    update_dict = GUPPY_MODEL_DESCRIPTIONS
    GUPPY_MODEL_DESCRIPTIONS = {important_models[0]:f"Default set by environmental variable"}
    GUPPY_MODEL_DESCRIPTIONS.update(update_dict)

CLAIR3_MODEL_DESCRIPTIONS = {
    "/opt/models/r941_prom_hac_g360+g422": "Default within clair3 Singularity container",
    "/opt/models/r941_prom_sup_g5014": "Super accurate model for v9.4.1 chemistry",
    "Other": "Choose your own adventure: specify below",
}

if os.getenv('CLAIR3_MODEL') is not None:
    update_dict = CLAIR3_MODEL_DESCRIPTIONS
    CLAIR3_MODEL_DESCRIPTIONS = {important_models[1]:f"Default set by environmental variable"}
    CLAIR3_MODEL_DESCRIPTIONS.update(update_dict)

max_threads = psutil.cpu_count(logical=True)
max_mem = round(bytesto(psutil.virtual_memory().available, "m"))


## DEFINE CLASSES ##
class AnalysisParameters:
    def __init__(self):
        self.barcode_kit = ""
        self.variant_caller = ""
        self.consensus_freq = 0.75
        self.indel_freq = 0.4
        self.reference = ""
        self.module = ""
        self.demultiplexed = ""
        # entry parameters
        self.samples_file = ""
        self.outdir = ""
        # advanced params
        self.snv_min_freq = ""
        self.min_depth = ""
        self.guppy_model = ""
        self.clair3_model = ""
        self.threads = ""
        self.max_memory = ""
        self.skip_clipping = ""
        # checkbox params
        self.demultiplexed = ""
        self.force = ""
        self.redo_analysis = ""
        self.delete_reads = ""
        self.print_dag = ""
        self.dry_run = ""
        self.create_envs_only = ""
        self.no_update = ""
        self.quiet = ""
        self.float_params = ['consensus_freq',
                             'indel_freq',
                             'snv_min_freq',
                             ]
        self.checkbox_params = ['demultiplexed',
                                'force',
                                'redo_analysis',
                                'delete_reads',
                                'print_dag',
                                'dry_run',
                                'create_envs_only',
                                'no_update',
                                'quiet',
                                ]  # Add the keys of the parameters derived from checkboxes

class AnalysisGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.parameters = AnalysisParameters()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("oat: ONT Analysis Toolkit")

        # main layout
        layout = QVBoxLayout()
        self.setLayout(layout)

        grid_layout = QGridLayout()
        layout.addLayout(grid_layout)

        # advanced options layout
        self.advanced_group = QGroupBox("Advanced Options")
        self.advanced_group.setCheckable(True)
        self.advanced_group.setChecked(False)
        layout.addWidget(self.advanced_group)
        advanced_layout = QGridLayout()
        self.advanced_group.setLayout(advanced_layout)

        ## MAIN OPTIONS ##
        #%% Barcode kit
        self.label_barcode_kit = QLabel("Barcode kit:")
        grid_layout.addWidget(self.label_barcode_kit, 0, 0)

        self.combo_barcode_kit = QComboBox()
        self.combo_barcode_kit.addItems(list(BARCODE_KIT_DESCRIPTIONS.keys()))
        self.combo_barcode_kit.setCurrentText("SQK-RBK004")
        self.combo_barcode_kit.currentTextChanged.connect(
            lambda option: self.update_desc(self.label_barcode_kit_desc, BARCODE_KIT_DESCRIPTIONS, option)
        )
        grid_layout.addWidget(self.combo_barcode_kit, 0, 1)

        self.label_barcode_kit_desc = QLabel("Rapid 12-barcode kit")
        grid_layout.addWidget(self.label_barcode_kit_desc, 0, 2)

        #%% Demultiplexing
        self.label_demultiplexed = QLabel("Demultiplexed:")
        grid_layout.addWidget(self.label_demultiplexed, 1, 0)

        # self.combo_demultiplexed = QComboBox()
        # self.combo_demultiplexed.addItems(["True", "False"])
        # self.combo_demultiplexed.setCurrentText("True")
        # self.combo_demultiplexed.currentTextChanged.connect(
        #     lambda option: self.update_desc(self.label_demultiplexed_desc, DEMULTIPLEXED_DESCRIPTIONS, option)
        # )
        # grid_layout.addWidget(self.combo_demultiplexed, 1, 1)

        # self.label_demultiplexed_desc = QLabel("Sample demultiplexing handled via MinKNOW")
        # grid_layout.addWidget(self.label_demultiplexed_desc, 1, 2)
        self.checkbox_demultiplexed = QCheckBox("")
        self.checkbox_demultiplexed.setChecked(True)
        grid_layout.addWidget(self.checkbox_demultiplexed, 1, 1)
        self.label_demultiplexed_desc = QLabel("Sample demultiplexing handled via MinKNOW")
        grid_layout.addWidget(self.label_demultiplexed_desc, 1, 2)

        # %% Reference
        self.label_reference = QLabel("Reference:")
        grid_layout.addWidget(self.label_reference, 2, 0)

        self.combo_reference = QComboBox()
        self.combo_reference.addItems(list(REFERENCE_DESCRIPTIONS.keys()))
        self.combo_reference.setCurrentText("MN908947.3")
        self.combo_reference.currentTextChanged.connect(
            lambda option: self.update_desc(self.label_reference_desc, REFERENCE_DESCRIPTIONS, option)
        )
        grid_layout.addWidget(self.combo_reference, 2, 1)

        self.label_reference_desc = QLabel("SARS-CoV-2 (Wuhan Hu-1)")
        grid_layout.addWidget(self.label_reference_desc, 2, 2)

        # %% Module
        self.label_module = QLabel("Module:")
        grid_layout.addWidget(self.label_module, 3, 0)

        self.combo_module = QComboBox()
        self.combo_module.addItems(["All", "Rampart", "Analysis"])
        self.combo_module.setCurrentText("True")
        self.combo_module.currentTextChanged.connect(
            lambda option: self.update_desc(self.label_module_desc, MODULE_DESCRIPTIONS, option)
        )
        grid_layout.addWidget(self.combo_module, 3, 1)

        self.label_module_desc = QLabel("Monitor run with RAMPART then run analysis")
        grid_layout.addWidget(self.label_module_desc, 3, 2)

        #%% Input spreadsheet widget
        self.label_samples_file = QLabel("Input Spreadsheet:")
        grid_layout.addWidget(self.label_samples_file, 4, 0)

        self.entry_samples_file = QLineEdit()
        grid_layout.addWidget(self.entry_samples_file, 4, 1)

        self.button_samples_file = QPushButton("Browse", self)
        self.button_samples_file.clicked.connect(
            lambda option: self.browse_file(self.entry_samples_file)
        )
        grid_layout.addWidget(self.button_samples_file, 4, 2)

        ## ADVANCED OPTIONS ##
        #%% snv min freq widget
        self.label_snv_min_freq = QLabel("SNV minimum frequency:")
        advanced_layout.addWidget(self.label_snv_min_freq, 0, 0)

        self.entry_snv_min_freq = QLineEdit("0.2")
        self.entry_snv_min_freq.textChanged.connect(partial(self.validate_float, "SNV minimum frequency"))
        advanced_layout.addWidget(self.entry_snv_min_freq, 0, 1)

        self.label_snv_min_freq_desc = QLabel("Minimum allele frequency for an SNV to be kept during variant calling.")
        advanced_layout.addWidget(self.label_snv_min_freq_desc, 0, 2)

        #%% Consensus SNP Frequency widget
        self.label_consensus_freq = QLabel("Consensus SNP frequency:")
        advanced_layout.addWidget(self.label_consensus_freq, 1, 0)

        self.entry_consensus_freq = QLineEdit("0")
        self.entry_consensus_freq.textChanged.connect(partial(self.validate_float, "Consensus SNP frequency"))

        advanced_layout.addWidget(self.entry_consensus_freq, 1, 1)

        self.label_consensus_freq_desc = QLabel("Variant allele frequency threshold for an SNP to be incorporated into consensus genome.\nVariants below this frequency will be incorporated with an IUPAC ambiguity.\nSet to 0 to incorporate the majority or most common base.\nNote: currently do not recommend anything except the default - debugging.\n")
        advanced_layout.addWidget(self.label_consensus_freq_desc, 1, 2)

        #%% Consensus Indel Frequency widget
        self.label_indel_freq = QLabel("Consensus indel frequency:")
        advanced_layout.addWidget(self.label_indel_freq, 2, 0)

        self.entry_indel_freq = QLineEdit("0.4")
        self.entry_indel_freq.textChanged.connect(partial(self.validate_float, "Consensus indel frequency"))
        advanced_layout.addWidget(self.entry_indel_freq, 2, 1)

        self.label_indel_freq_desc = QLabel("Variant allele frequency threshold for an indel to be incorporated into consensus genome.\nVariants below this frequency will not be incorporated.\nSet to 0 to incorporate the majority or most common base.")
        advanced_layout.addWidget(self.label_indel_freq_desc, 2, 2)

        #%% Minimum depth
        self.label_min_depth = QLabel("Minimum depth:")
        advanced_layout.addWidget(self.label_min_depth, 3, 0)

        self.entry_min_depth = QLineEdit("15")
        advanced_layout.addWidget(self.entry_min_depth, 3, 1)

        self.label_min_depth_desc = QLabel("Minimum depth for (1) an SNV to be kept and (2) consensus genome generation.")
        advanced_layout.addWidget(self.label_min_depth_desc, 3, 2)

        #%% variant caller
        self.label_variant_caller = QLabel("Variant caller:")
        advanced_layout.addWidget(self.label_variant_caller, 4, 0)

        self.combo_variant_caller = QComboBox()
        self.combo_variant_caller.addItems(['clair3', 'medaka'])
        self.combo_variant_caller.setCurrentText("clair3")
        advanced_layout.addWidget(self.combo_variant_caller, 4, 1)

        #%% guppy model
        self.label_guppy_model = QLabel("Guppy model:")
        advanced_layout.addWidget(self.label_guppy_model, 5, 0)

        self.combo_guppy_model = QComboBox()
        self.combo_guppy_model.addItems(list(GUPPY_MODEL_DESCRIPTIONS.keys()))
        self.combo_guppy_model.setCurrentText(str(important_models[0]))
        advanced_layout.addWidget(self.combo_guppy_model, 5, 1)
        self.label_guppy_model_desc = QLabel("Select a guppy model or select 'Other' and specify")
        self.combo_guppy_model.currentTextChanged.connect(
            lambda option, label=self.label_guppy_model_desc, input_dict=GUPPY_MODEL_DESCRIPTIONS: self.update_desc(label, input_dict, option)
        )
        advanced_layout.addWidget(self.label_guppy_model_desc, 5, 2)

        self.entry_other = QLineEdit()
        advanced_layout.addWidget(self.entry_other, 6, 1)

        self.combo_guppy_model.currentTextChanged.connect(self.update_other_entry)

        #%% clair3 model
        self.label_clair3_model = QLabel("Clair3 model:")
        advanced_layout.addWidget(self.label_clair3_model, 7, 0)

        self.combo_clair3_model = QComboBox()
        self.combo_clair3_model.addItems(list(CLAIR3_MODEL_DESCRIPTIONS.keys()))
        self.combo_clair3_model.setCurrentText(str(important_models[1]))
        advanced_layout.addWidget(self.combo_clair3_model, 7, 1)
        self.label_clair3_model_desc = QLabel("Select a clair3 model or select 'Other' and specify")
        self.combo_clair3_model.currentTextChanged.connect(
            lambda option, label=self.label_clair3_model_desc, input_dict=CLAIR3_MODEL_DESCRIPTIONS: self.update_desc(label, input_dict, option)
        )
        advanced_layout.addWidget(self.label_clair3_model_desc, 7, 2)

        self.entry_other = QLineEdit()
        advanced_layout.addWidget(self.entry_other, 8, 1)

        self.combo_clair3_model.currentTextChanged.connect(self.update_other_entry)

        #%% threads
        self.label_threads = QLabel("Number of threads:")
        advanced_layout.addWidget(self.label_threads, 9, 0)

        self.entry_threads = QLineEdit(str(max_threads))
        advanced_layout.addWidget(self.entry_threads, 9, 1)

        #%% max memory
        self.label_max_memory = QLabel("Maximum memory (in MB):")
        advanced_layout.addWidget(self.label_max_memory, 10, 0)

        self.entry_max_memory = QLineEdit(str(max_mem))
        advanced_layout.addWidget(self.entry_max_memory, 10, 1)

        self.label_max_memory_desc = QLabel("Defaults to using most of your available RAM")
        advanced_layout.addWidget(self.label_max_memory_desc, 10, 2)

        #%% ADD OTHER OPTIONS HEADING
        self.label_other_options = QLabel("Other options:")
        advanced_layout.addWidget(self.label_other_options, 11, 0)

        #%% overwrite analysis
        self.checkbox_force = QCheckBox("Overwrite previous analysis from Snakemake stage")
        advanced_layout.addWidget(self.checkbox_force, 12, 1)

        #%% redo analysis
        self.checkbox_redo_analysis = QCheckBox("Redo previous analysis from scratch")
        advanced_layout.addWidget(self.checkbox_redo_analysis, 13, 1)

        #%% delete reads
        self.checkbox_delete_reads = QCheckBox("Delete demultiplexed reads after analysis")
        advanced_layout.addWidget(self.checkbox_delete_reads, 14, 1)

        #%% skip clipping
        self.label_skip_clipping = QLabel("Skip primer clipping:")
        advanced_layout.addWidget(self.label_skip_clipping, 15, 0)

        self.combo_skip_clipping = QComboBox()
        self.combo_skip_clipping.addItems(list(SKIP_CLIPPING_DESCRIPTIONS.keys()))
        self.combo_skip_clipping.setCurrentText(False)
        advanced_layout.addWidget(self.combo_skip_clipping, 15, 1)
        self.label_skip_clipping_desc = QLabel("Choose whether amplicon primer skipping should be skipped")
        self.combo_skip_clipping.currentTextChanged.connect(
            lambda option, label=self.label_skip_clipping_desc, input_dict=SKIP_CLIPPING_DESCRIPTIONS: self.update_desc(label, input_dict, option)
        )
        advanced_layout.addWidget(self.label_skip_clipping_desc, 15, 2)

        self.entry_other = QLineEdit()
        advanced_layout.addWidget(self.entry_other, 15, 1)

        self.combo_skip_clipping.currentTextChanged.connect(self.update_other_entry)

        #%% print DAG
        self.checkbox_print_dag = QCheckBox("Print analysis DAG then quit")
        advanced_layout.addWidget(self.checkbox_print_dag, 15, 1)

        #%% dry run
        self.checkbox_dry_run = QCheckBox("Dry run only")
        advanced_layout.addWidget(self.checkbox_dry_run, 16, 1)

        #%% create envs only
        self.checkbox_create_envs_only = QCheckBox("Create conda envs only")
        advanced_layout.addWidget(self.checkbox_create_envs_only, 17, 1)

        #%% no update
        self.checkbox_no_update = QCheckBox("Disable updating of container versions for SARS-CoV-2 analysis.")
        advanced_layout.addWidget(self.checkbox_no_update, 18, 1)

        #%% quiet
        self.checkbox_quiet = QCheckBox("Stop printing of Snakemake commands to screen")
        advanced_layout.addWidget(self.checkbox_quiet, 19, 1)

        #%% Outdir
        self.label_outdir = QLabel("Output directory:")
        advanced_layout.addWidget(self.label_outdir, 20, 0)

        self.entry_outdir = QLineEdit()
        advanced_layout.addWidget(self.entry_outdir, 20, 1)

        self.button_outdir = QPushButton("Browse", self)
        self.button_outdir.clicked.connect(
            lambda option: self.browse_directory(self.entry_outdir)
        )
        advanced_layout.addWidget(self.button_outdir, 20, 2)

        ## RUN ANALYSIS ##
        #%% Run analysis button
        self.button_run_analysis = QPushButton("Run Analysis", self)
        self.button_run_analysis.setStyleSheet("background-color: lightpink; font-weight: bold;")
        self.button_run_analysis.clicked.connect(self.run_analysis)
        layout.addWidget(self.button_run_analysis)

    def update_other_entry(self, selected_option):
        if selected_option == 'Other':
            self.entry_other.setEnabled(True)
        else:
            self.entry_other.setEnabled(False)
            self.entry_other.clear()

    def update_desc(self, label_desc, input_dict, option):
        description = input_dict[option]
        label_desc.setText(description)

    def browse_file(self, entry):
        filename, _ = QFileDialog.getOpenFileName()
        entry.setText(filename)

    def browse_directory(self, entry):
        directory = QFileDialog.getExistingDirectory()
        entry.setText(directory)

    def check_parameters(self):
        if not 0 <= float(self.entry_consensus_freq.text()) <= 1:
            QMessageBox.warning(self, "Invalid consensus SNP frequency",
                                "Value must be between 0 and 1.")
            return False
        if not 0 <= float(self.entry_indel_freq.text()) <= 1:
            QMessageBox.warning(self, "Invalid consensus indel frequency",
                                "Value must be between 0 and 1.")
            return False
        if not 0 <= float(self.entry_snv_min_freq.text()) <= 1:
            QMessageBox.warning(self, "Invalid SNV minimum frequency",
                                "Value must be between 0 and 1.")
            return False
        if int(self.entry_threads.text()) < 1:
            QMessageBox.warning(self, "Invalid number of threads",
                                "Value must be >= 1.")
            return False
        if max_threads < int(self.entry_threads.text()):
            QMessageBox.warning(self, "Too many threads",
                                f"Value should not exceed the maximum threads available ({max_threads}).")
            return False
        if not self.entry_samples_file.text():
            QMessageBox.critical(self, "Error", "Hmm, forgetting something? You need to choose an input spreadsheet.")
            return False
        QMessageBox.information(self, "Run successfully deployed",
                                "Run deployed. Check the Terminal to monitor the analysis.")
        return True

    def generate_cli_input(self):
        cli_input = []
        for key, value in vars(self.parameters).items():
            if key in ["checkbox_params","float_params"]:
                continue
            elif key in self.parameters.checkbox_params:
                if value:
                    cli_input.append('--' + str(key))
            elif key != "samples_file":
                cli_input.append('--' + str(key))
                cli_input.append(str(value))
            elif key == "samples_file":
                cli_input.append(value)
            else:
                print("Something went wrong")
        return cli_input

    def run_analysis(self):
        self.update_parameters()
        success = self.check_parameters()
        if success:
            self.close()
            # Add your code here to run the analysis with the updated parameters
            import oat.cli.cli as cli
            import sys
            # generate cli input
            sys.argv = [sys.argv[0]] + [str(x) for x in self.generate_cli_input()]
            cli.main(sys.argv[1:])

    def validate_float(self, title, text):
        try:
            value = float(text)
            if not 0 <= value <= 1:
                QMessageBox.warning(self, f"Invalid {title}", "Value must be between 0 and 1.")
        except ValueError:
                pass
    def update_parameters(self):
        # main options
        self.parameters.barcode_kit = self.combo_barcode_kit.currentText()
        self.parameters.reference = self.combo_reference.currentText()
        self.parameters.module = self.combo_module.currentText()
        # advanced options
        self.parameters.variant_caller = self.combo_variant_caller.currentText()
        self.parameters.consensus_freq = float(self.entry_consensus_freq.text())
        self.parameters.snv_min_freq = float(self.entry_snv_min_freq.text())
        self.parameters.indel_freq = float(self.entry_indel_freq.text())
        self.parameters.min_depth = int(self.entry_min_depth.text())
        self.parameters.guppy_model = self.combo_guppy_model.currentText()
        self.parameters.clair3_model = self.combo_clair3_model.currentText()
        self.parameters.threads = int(self.entry_threads.text())
        self.parameters.max_memory = int(self.entry_max_memory.text())
        self.parameters.outdir = self.entry_outdir.text()

        # checkbox parameters = flag-only options ('store_true')
        self.parameters.demultiplexed = self.checkbox_demultiplexed.isChecked()
        self.parameters.force = self.checkbox_force.isChecked()
        self.parameters.redo_analysis = self.checkbox_redo_analysis.isChecked()
        self.parameters.delete_reads = self.checkbox_delete_reads.isChecked()
        self.parameters.print_dag = self.checkbox_print_dag.isChecked()
        self.parameters.dry_run = self.checkbox_dry_run.isChecked()
        self.parameters.create_envs_only = self.checkbox_create_envs_only.isChecked()
        self.parameters.no_update = self.checkbox_no_update.isChecked()
        self.parameters.quiet = self.checkbox_quiet.isChecked()
        # positional argument
        self.parameters.samples_file = self.entry_samples_file.text()


def main():
    app = QApplication(sys.argv)
    analysis_gui = AnalysisGUI()
    analysis_gui.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
