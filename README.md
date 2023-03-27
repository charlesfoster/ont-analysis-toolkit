[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
```
                                       ,d
                                       88
              ,adPPYba,  ,adPPYYba, MM88MMM
              a8"     "8a ""     `Y8   88
              8b       d8 ,adPPPPP88   88
              "8a,   ,a8" 88,    ,88   88,
              `"YbbdP"'  `"8bbdP"Y8   "Y888

```
# ONT Analysis Toolkit (OAT)
A pipeline to facilitate sequencing of viral genomes (amplified with a tiled amplicon scheme) and assembly into consensus genomes. Supported viruses currently include Human betaherpesvirus 5 (CMV) and SARS-CoV-2. ONT sequencing data from a MinION can be monitored in real time with the `rampart` module. All analysis steps are handled by the `analysis` module, whereby sequencing data are analysed using a pipeline written in snakemake, with the choice of tools heavily influenced by the `artic minion` pipeline. Both steps can be run in order with a single command using the `all` module.

Author: Dr Charles Foster

# Starting out
To begin with, clone this github repository:

```
git clone https://github.com/charlesfoster/ont-analysis-toolkit.git

cd ont-analysis-toolkit
```

Next, install *most* dependencies using `conda`:

```
conda env create -f environment.yml
```

Pro tip: if you install `mamba`, you can create the environment with that command instead of `conda`. A lot of `conda` headaches go away: it's a much faster drop-in replacement for `conda`.

```
conda install mamba
mamba env create -f environment.yml
```

Install the pipeline correctly by activating the conda environment and using the `setup.py` script with:

```
conda activate oat
pip install .
```

Other dependencies:
* Demultiplexing is done using `guppy_barcoder`. The program will need to be installed and in your path.
* If analysing SARS-CoV-2 data, lineages are typed using `pangolin`. Accordingly, `pangolin` needs to be installed according to instructions at https://github.com/cov-lineages/pangolin. The pipeline will fail if the `pangolin` environment can't be activated.
* Variants are called using `medaka` and `longshot`. Ideally we could install these via `mamba` in the main `environment.yml` file, but there are sadly some necessary libraries for `medaka` that are incompatible with our main `oat` environment. Consequently, I have written the analysis pipeline so that `snakemake` automatically installs `medaka` and its dependencies into an isolated environment during execution of the `oat` pipeline. The environment is only created the first time you run the pipeline, *but* if you run the pipeline from a different working directory in the future, the environment will be created again. Solution: always run `oat` from the same working directory

**tl;dr**: you don't need to do anything for variant calling to work; just don't get confused during the initial pipeline run when the terminal indicates creation of a new conda environment. You should run `oat` from the same directory each time, otherwise a new conda environment will be created each time.

# Usage
The environment with all necessary tools is installed as '`oat`' for brevity. The environment should first be activated:

```
conda activate oat
```

Then, to run the pipeline, it's as simple as:

```
oat <input_spreadsheet.csv>
```

where <input_spreadsheet.csv> should be replaced with the full path to a spreadsheet with minimal metadata for the ONT sequencing run (see example spreadsheet: run_data_example.csv). The most important thing to remember is that the 'run_name' in the spreadsheet must exactly match the name of the sequencing run, as set up in MinKNOW.

Note that there are many additional options/settings to take advantage of:

```

                                        ,d
                                        88
               ,adPPYba,   ,adPPYYba, MM88MMM
              a8"     "8a ""     `Y8    88
              8b       d8 ,adPPPPP88    88
              "8a,   ,a8" 88,    ,88    88,
               `"YbbdP"'  `"8bbdP"Y8   "Y888

        OAT: ONT Analysis Toolkit (version 0.10.5)

usage: oat [options] <samples_file>

A pipeline for sequencing and analysis of viral genomes using an ONT MinION

positional arguments:
  samples_file          Path to file with sample metadata (.csv format). See example spreadsheet for minimum necessary
                        information.

optional arguments:
  -h, --help            show this help message and exit
  -b BARCODE_KIT, --barcode_kit BARCODE_KIT

                                    Barcode kit that you used: necessary for demultiplexing.

                                    Common options include:
                                    - rapid 12-barcode kit (SQK-RBK004)
                                    - rapid 96-barcode kit (SQK-RBK110-96)
                                    - native 12-barcode kit (EXP-NBD104)
                                    - native 12-barcode expansion kit 13-24 (EXP-NBD114)
                                    - native 96-barcode kit (EXP-NBD196)

                                    Please ensure your version of guppy_barcoder supports the kit name. Try:
                                    'guppy_barcoder --print_kits'

                                    Default: SQK-RBK004

  -c <float>, --consensus_freq <float>

                                    Variant allele frequency threshold for a variant to be incorporated into consensus
                                    genome.
                                    Variants below this frequency will be incorporated with an IUPAC ambiguity.
                                    Set to 0 to incorporate the majority or most common base.
                                    Note: currently do not recommend anything except the default - debugging.
                                    Default: 0.0

  -i <float>, --indel_freq <float>

                                    Variant allele frequency threshold for an indel variant to be incorporated into
                                    consensus genome.
                                    Variants below this frequency will not be incorporated.
                                    Set to 0 to incorporate the majority or most common base.
                                    Default: 0.4

  -d, --demultiplexed    
                                    Reads already demultiplexed using guppy_barcoder into
                                    '/var/lib/minknow/data/<run_name>'.
                                    By default, assumes reads need to be demultiplexed and reads are demultiplexed
                                    into the output directory.

  -f, --force           Force overwriting of completed files in snakemake analysis (Default: files not overwritten)
  -n, --dry_run         Dry run only
  -m rampart | analysis | all, --module rampart | analysis | all
                        Pipeline module to run: 'rampart', 'analysis' or 'all' (rampart followed by analysis).
                        Default: 'all'
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: /home/vrl/Programs/ont-analysis-toolkit/analysis_results +
                        'organism_name' + 'run_name' from samples spreadsheet
  --rampart_outdir RAMPART_OUTDIR
                        Output directory. Default: /home/vrl/Programs/ont-analysis-toolkit/rampart_files
  -p, --print_dag       Save directed acyclic graph (DAG) of workflow to outdir
  -r REFERENCE, --reference REFERENCE

                                    Reference genome to use. Supported references:
                                        - 'MN908947.3' (SARS-CoV-2)
                                        - 'NC_006273.2' (CMV Merlin) (planned support, but not yet implemented)
                                    Other references can be used, but the corresponding assembly (fasta) and
                                    annotation (gff3 from Ensembl) must be added. See README.
                                    Default: MN908947.3

  -t <int>, --threads <int>
                        Number of threads to use
  -v VARIANT_CALLER, --variant_caller VARIANT_CALLER
                        Variant caller to use. Choices: 'clair3','medaka'. Default: 'clair3'
  --create_envs_only    Create conda environments for snakemake analysis, but do no further analysis. Useful for
                        initial pipeline setup. Default: False
  --snv_min_freq SNV_MIN_FREQ
                        Minimum variant allele frequency for an SNV to be kept. Default: 0.2
  --guppy_model GUPPY_MODEL
                        Model used within guppy for basecalling - needed for medaka analyses. 'GUPPY_MODEL'
                        environmental model is used, if set, otherwise default: r941_min_high_g360
  --min_depth MIN_DEPTH
                        Minimum depth for (1) an SNV to be kept; and (2) consensus genome generation. Default: 20
  --alternate_analysis  Run an alternate analysis to generate consensus genomes based on different input parameters.
                        Use after initial run.
  --alt_cov_max ALT_COV_MAX
                        Samples with a genome coverage alt_cov_min <= x < alt_cov_max will be chosen for an alternate
                        analysis. ONLY USED WITH ALTERNATE ANALYSIS OPTION. Default: 80
  --alt_cov_min ALT_COV_MIN
                        Samples with a genome coverage alt_cov_min <= x < alt_cov_max will be chosen for an alternate
                        analysis. ONLY USED WITH ALTERNATE ANALYSIS OPTION. Default: 40
  --delete_reads        Delete demultiplexed reads after analysis
  --redo_analysis       Delete entire analysis output directory and contents for a fresh run
  --version             show program's version number and exit
  --minknow_data MINKNOW_DATA
                        Location of MinKNOW data root. Default: /var/lib/minknow/data
  --max_memory <int>    Maximum memory (in MB) that you would like to provide to snakemake. Default: 53616MB
  --quiet               Stop printing of snakemake commands to screen.
  --report              Generate report (currently minimally functional).

Amplicon protocol must be specified in input spreadsheet.
```

# What does the pipeline do?
- [Monitoring of sequencing with RAMPART](#RAMPART-Module)
- [Analysis of sequencing data](#Analysis-Module)
- [Other notes](#Other-Notes)
- [Credits](#Credits)

# RAMPART Module
All input files for RAMPART are generated based on your input spreadsheet, and a web browser is launched to view the sequencing in real time.

# Analysis Module
Reads are mapped to the relevant reference genome with `minimap2`. Amplicon primers are trimmed using `samtools ampliconclip`. Variants are called using either `clair3` or `medaka` and `longshot`, followed by filtering and consensus genome assembly using `bcftools`. The amino acid consequences of SNPs are inferred using `bcftools csq`. If analysing SARS-CoV-2, lineages are inferred using `pangolin` and `nextclade`. Finally, a variety of sample QC metrics are combined into a final QC file.

# Other Notes
**Run metadata spreadsheet**
The run metadata spreadsheet should be styled after provided csv file ('run_data_example.csv'), with the minimum columns to include "barcode", "id", "neg_control", "protocol", and "run_name". Regardless of whether you use the ONT "native" (ligation) barcodes or the "rapid" (transposase-based) barcodes, the barcodes in the spreadsheet should be named in the form of "BCXX", where "XX" represents a number from 01 onwards.

**Barcodes**
By default, the pipeline assumes you are using the ONT 12-barcode rapid kit ("SQK-RBK004"). However, other kits are supported, including:
- SQK-LSK109 (up to 12 samples ligation/native barcodes)
- SQK-RBK110-96 (up to 96 samples rapid barcodes)

**Protocols**
A protocol for the amplicon scheme needs (a) to be installed in the pipeline, and (b) named in the run_data.csv spreadsheet for analyses to work correctly. The pipeline comes with the Midnight protocol for SARS-CoV-2 pre-installed (https://www.protocols.io/view/sars-cov2-genome-sequencing-protocol-1200bp-amplic-bwyppfvn). Adding additional protocols is fairly easy:

1. Make a directory called /path/to/ont-analysis-toolkit/oat/protocols/ARTICV3 (needs to be in all caps)
2. Make a directory within ARTICV3 called 'rampart'

   (a) Put the normal rampart files within that directory (genome.json, primers.json, protocol.json, references.fasta)
3. Make a directory within ARTICV3 called 'schemes'

    (a) Put the 'scheme.bed' file with primer coordinates in the 'schemes' directory

5. Edit the minimum and maximum length parameters in /path/to/ont-analysis-toolkit/oat/protocols/length_params.csv

6. Make sure you're in /path/to/ont-analysis-toolkit/, then activate the conda environment and use the following command: `pip install .`

Done!

Note: I've provided two script to help generate the `genome.json` file:
- /path/to/ont-analysis-toolkit/oat/scripts/genbank2json.py: converts a genbank-format file into the necessary json format for RAMPART
- /path/to/ont-analysis-toolkit/oat/scripts/gff2json.py: converts a gff3-format file into the necessary json format for RAMPART. Works from limited testing with either NCBI or Ensembl gff3 files.

**Amino acid consequences**
For the amino acid consequences step to work, a requirement is an annotation file for the chosen reference genome. The annotations must be in gff3 format, and must be in the 'Ensembl flavour' of gff3. There is a script included in the repository that can convert an NCBI gff3 file into an 'Ensembl flavour' gff3 file: `/path/to/ont-analysis-toolkit/oat/scripts/gff2gff.py`.

**Adding another reference organism**
To add a reference organism/assembly other than the provided references for SARS-CoV-2 and CMV Merlin, a reference genome is needed (fasta format), as well as an annotation file in 'Ensembl flavour' gff3 format (see [section on amino acid consequences](#Amino-acid-consequences)). Both files need to be added to `/path/to/ont-analysis-toolkit/oat/references`, and should be named consistently, e.g. `my_reference.fasta` and `my_reference.gff3`. After adding the files, make sure you're in /path/to/ont-analysis-toolkit/, then activate the conda environment and use the following command: `pip install .` The next time you run the pipeline, you can then specify your new reference like so: `oat -r new_reference input_spreadsheet.csv`.

# Credits
* When this pipeline is used, citations should be found for the programs used internally.
* The gff3 file I included for SARS-CoV-2 was originally sent to me by Torsten Seemann.
* Being new to using snakemake + wrapper scripts, I used `pangolin` as a guide for directory structure and rule creation - so thanks to them.
* The analysis module was heavily influenced by the ARTIC team, especially the `artic minion` pipeline.
* `gb2gff.py` is based on work by Damien Farrell https://dmnfarrell.github.io/bioinformatics/bcftools-csq-gff-format
