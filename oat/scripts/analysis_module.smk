#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last edited on Thur December 02 2021

@author: Dr Charles Foster
"""

###########
# Libraries
###########

import glob
import re
from datetime import date
import pandas as pd
import yaml
import subprocess
import glob
import pathlib
import pandas as pd
import os

#################
# Custom functions
#################


###############
# Configuration
###############
# set up env variable & threads for tensorflow / medaka
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
#medaka_con_threads = int(config["threads"] / 2)
# medaka log consistently says more than 2 threads is a waste, so...
medaka_con_threads = int(2)

# define some key variables
RESULT_DIR = config["outdir"]
SAMPLES = [
    os.path.basename(x).replace(".fastq", "")
    for x in glob.glob(config["reads_dir"] + "/*.fastq")
]

TODAY = date.today().strftime("%Y-%m-%d")
if os.path.basename(config["reference"]) == "MN908947.3.fasta":
    SARS_ANALYSIS = True
else:
    SARS_ANALYSIS = False

# index reference genomes if necessary
if not os.path.isfile(config["reference"] + ".fai"):
    os.system("samtools faidx {} 2> /dev/null".format(config["reference"]))
if not os.path.isfile(config["reference"] + ".bwt"):
    os.system("bwa index {} 2> /dev/null".format(config["reference"]))

# conda env for lofreq - not used currently
if config["variant_caller"] == "lofreq":
    variant_conda_env = "../envs/lofreq.yaml"

# set up column name for coverage - dependent on input parameter
# update: for simplifying reporting purposes, changing to just 'coverage'
# depth value will still be in the config
config["min_depth"] = int(config["min_depth"])
#coverage_colname = "ref_cov_"+str(config["min_depth"])
coverage_colname = "coverage"

# determine model for clair3
base_gmodel = config["guppy_model"]

if config["variant_caller"].lower() == "clair3":
    if '_g5' in base_gmodel:
        clair3_model = 'r941_prom_sup_g5014'
    else:
        print('Could not determine appropriate clair3 model - setting to r941_prom_hac_g360+g422 for safety')
        clair3_model = 'r941_prom_hac_g360+g422'
else:
    clair3_model = 'None'

if config['variant_caller'] == 'clair3':
    snv_min_qual = 5
    filter_extension = "clair3.vcf.gz"
elif config['variant_caller'] == 'medaka':
    snv_min_qual = 20
    filter_extension = "longshot.vcf.gz"
elif config['variant_caller'] == 'lofreq':
    snv_min_qual = 20
    filter_extension = "lofreq.vcf.gz"

# medaka script
snakedir = os.path.abspath(os.path.dirname(__file__))
one_up = os.path.split(snakedir)[0]

parse_script = os.path.abspath(os.path.join(one_up,"oat","scripts","parse_medaka_variants.py"))

################
# Optional removal of trimmed reads (to save space)
################


onsuccess:
    if config["delete_reads"] == True:
        print(
            "\n\033[92mRemoving trimmed reads (input reads remain untouched)\033[0m\n"
        )
        dead_reads = glob.glob(config["reads_dir"] + "*.fastq", recursive=False)
        [os.remove(x) for x in dead_reads]

    print("\n\033[92mRemoving unwanted and/or empty files\033[0m\n")
    for file in glob.glob(RESULT_DIR + "/**", recursive=True):
        if file.endswith(
            (
                "draft.vcf.gz.tbi",
                ".filtered.vcf.gz.csi",
                ".mapped.bam.csi",
            )
        ) or (os.path.getsize(file) == 0 and not file.endswith(".fastq")):
            os.remove(file)


################
# rules
################

potential_pangolin = []
potential_nextclade = []
if SARS_ANALYSIS:
    potential_pangolin.append(
        expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.lineage_report.csv"),
            sample=SAMPLES,
        )
    )
    potential_nextclade.append(
        expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.nextclade_report.tsv"),
            sample=SAMPLES,
        )
    )
else:
    potential_pangolin.append(
        expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.qc_results.csv"),
            sample=SAMPLES,
        )
    )
    potential_nextclade.append(
        expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
            sample=SAMPLES,
        )
    )

rule final_qc:
    input:
        consensus_genome=expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.consensus.fasta"),
            sample=SAMPLES,
        ),
        pangolin_results=potential_pangolin,
        nextclade_results=potential_nextclade,
        bcsq_tsv=expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.annotated.tsv"), sample=SAMPLES
        ),
        reports=expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.qc_results.csv"),
            sample=SAMPLES,
        ),
    params:
        run_metadata=os.path.join(RESULT_DIR, "metadata.csv"),
        sars_analysis=SARS_ANALYSIS,
    run:
        # combine fasta files into one multifasta
        fa_files = [f for f in glob.glob(RESULT_DIR + "/**/*.consensus.fasta", recursive=True) if (s in f for s in SAMPLES)]

        multifasta = os.path.join(
            RESULT_DIR, config["run_name"] + ".consensus_genomes.fa"
        )
        if os.path.exists(multifasta):
            os.remove(multifasta)

        for fa in fa_files:
            os.system("cat {0} >> {1}".format(fa, multifasta))
        run_metadata = pd.read_csv(params.run_metadata)

        # collect QC files
        qc_files = []
        for sample in SAMPLES:
            sample_dir = os.path.join(RESULT_DIR, sample)
            qc_file = glob.glob(sample_dir + "/*qc_results.csv")
            if len(qc_file) == 1:
                qc_files.append(qc_file[0])
            else:
                print(f'Error: investigate duplicated sample_qc file found for {sample}')
                sys.exit(-1)

        # combine into one final_qc
        #qc_files = glob.glob(RESULT_DIR + "/**/*qc_results.csv", recursive=True)
        combined_qc = pd.concat([pd.read_csv(f) for f in qc_files]).set_index(["id"])
        outdata = run_metadata.join(combined_qc, on=["id"])
        outdata = outdata.set_index("id", drop=False)
        outdata.index.name = None

        # add analysis date
        outdata["analysis_date"] = date.today().strftime("%Y-%m-%d")

        # define order of columns in final_qc file
        if params.sars_analysis:

            # read in the lineage files
            lineages = pd.concat([pd.read_csv(f) for f in input.pangolin_results]).drop(
                ["ambiguity_score", "scorpio_conflict", "scorpio_notes", "is_designated", "qc_notes"], axis=1
            )

            lineages.columns = [
                "id",
                "lineage",
                "lineage_conflict",
                "scorpio_call",
                "scorpio_support",
                "lineage_designation_version",
                "pangolin_version",
                "scorpio_version",
                "constellation_version",
                "pangolin_status",
                "pangolin_note",
            ]

            outdata = outdata.join(lineages.set_index(["id"]), on=["id"])

            # read in the nextclade files
            nextclades = pd.concat([pd.read_csv(f, sep="\t") for f in input.nextclade_results])
            nextclades.rename(columns={'seqName':'id','qc.privateMutations.total':'totalPrivateMutations','qc.overallStatus':'Nextclade_QC'}, inplace=True)
            keep = ['id','clade','Nextclade_pango','Nextclade_QC','totalFrameShifts','totalAminoacidInsertions','totalAminoacidDeletions','totalAminoacidSubstitutions','totalNonACGTNs','totalPrivateMutations']
            nextclades = nextclades.loc[:,keep]

            outdata = outdata.join(nextclades.set_index(["id"]), on=["id"])

            compulsory_col_order = [
                "analysis_date",
                "run_name",
                "id",
                "barcode",
                "lineage",
                "scorpio_call",
                "Nextclade_pango",
                "clade",
                coverage_colname,
                "mean_depth",
                "num_reads",
                "num_mapped_reads",
                "percent_total_reads",
                "reads_qc",
                "pangolin_status",
            ]
        else:
            compulsory_col_order = [
                "analysis_date",
                "run_name",
                "id",
                "barcode",
                coverage_colname,
                "mean_depth",
                "num_reads",
                "num_mapped_reads",
                "percent_total_reads",
                "reads_qc",
            ]

        # do qc of number of reads
        reads_qc = pd.DataFrame(columns=["percent_total_reads", "reads_qc"])
        outdata = outdata.join(reads_qc, how="outer")

        total_reads = combined_qc["num_reads"].sum()
        sample_dict = dict(tuple(outdata.groupby("id")))
        for sample in sample_dict:
            percent = (sample_dict[sample]["num_reads"].squeeze() / total_reads) * 100
            outdata.at[sample, "percent_total_reads"] = percent
            if sample_dict[sample]["neg_control"].bool() == True:
                qc_result = "PASS" if percent < 5 else "FAIL"
            else:
                qc_result = "PASS" if percent > 5 else "FAIL"
            outdata.at[sample, "reads_qc"] = qc_result

        input_cols = outdata.columns.to_list()
        extra_input_cols = list(set(input_cols) - set(compulsory_col_order))
        final_col_order = compulsory_col_order + extra_input_cols
        outdata = outdata[final_col_order]
        outdata.to_csv(
            os.path.join(RESULT_DIR, config["run_name"] + "_qc.csv"),
            index=False,
            na_rep="NA",
        )
        os.remove(params.run_metadata)
        with open(os.path.join(RESULT_DIR, "config.yaml"), "w") as outfile:
            yaml.dump(config, outfile, default_flow_style=False)



rule map_reads:
    input:
        fastq=os.path.join(config["reads_dir"], "{sample}.fastq"),
    output:
        #bam=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.mapped.bam")),
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.mapped.bam"),
    message:
        "mapping reads for {wildcards.sample} to reference"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/minimap2.log.txt"),
    params:
        reference=config["reference"],
    resources:
        cpus=4,
    shell:
        """
        minimap2 -a -x map-ont -t {threads} {params.reference} {input.fastq} 2>{log} | \
        samtools view -@{threads} -bS -F 4 - | \
        samtools sort -@{threads} --write-index -o {output.bam} - 2>>{log}
        """


rule trim_amplicon_primers:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.mapped.bam"),
    output:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    message:
        "trimming amplicon primers from {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/ampliconclip.log.txt"),
    params:
        reference=config["reference"],
        bed=config["protocol_bed"],
    resources:
        cpus=4,
    shell:
        """
        samtools ampliconclip --both-ends --strand --soft-clip --filter-len 30 -@ 4 -u -b {params.bed} {input.bam} 2>{log} | \
        samtools sort -u -@ 4 -n 2>>{log} | \
        samtools fixmate -u -@ 4 - - 2>>{log} | \
        samtools calmd -u -@ 4 - {params.reference} 2>>{log} | \
        samtools sort -u -@ 4 2>>{log} | \
        samtools view --write-index -@ 20 -F 4 -o {output.bam}
        """

##### MEDAKA VARIANT CALLING START #####
rule medaka_consensus:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        hdf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.hdf")),
    message:
        "initial medaka consensus for {wildcards.sample}"
    threads: medaka_con_threads
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/medaka_consensus.log.txt"),
    params:
        reference=config["reference"],
        model=config["guppy_model"],
    resources:
        cpus=medaka_con_threads,
        gpu=1,
    conda:
        "../envs/medaka.yaml"
    shell:
        """
        export TF_FORCE_GPU_ALLOW_GROWTH=true; medaka consensus --model {params.model} --threads {threads} --chunk_len 800 --chunk_ovlp 400 {input.bam} {output.hdf} 2>{log}
        """


rule medaka_variant:
    input:
        hdf=os.path.join(RESULT_DIR, "{sample}/{sample}.hdf"),
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        medaka_vcf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf")),
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf.gz"),
    message:
        "initial medaka variant calls for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/medaka_variant.log.txt"),
    params:
        reference=config["reference"],
        model=config["guppy_model"],
    resources:
        cpus=4,
        gpu=1,
    conda:
        "../envs/medaka.yaml"
    shell:
        """
        export TF_FORCE_GPU_ALLOW_GROWTH=true; medaka variant {params.reference} {input.hdf} {output.medaka_vcf} 2>{log}
        bgzip -c {output.medaka_vcf} > {output.vcf} 2>>{log}
        tabix -p vcf {output.vcf} 2>>{log}
        """
        #medaka tools annotate --pad 1 --chunk_size 29903 --dpsp {output.medaka_vcf} {params.reference} {input.bam} {output.vcf} 2>>{log}

rule longshot:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf.gz"),
    output:
        #vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.all.vcf"),
        vcf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.all.vcf")),
    message:
        "neater longshot variant calls for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/longshot.log.txt"),
    params:
        reference=config["reference"],
        model="r941_min_high_g360",
    resources:
        cpus=4,
    shell:
        """
         longshot -P 0 -F --no_haps --bam {input.bam} --ref {params.reference} --out {output.vcf} --potential_variants {input.vcf}    2>{log}
        """


rule add_rough_VAF:
    input:
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.all.vcf"),
    output:
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.longshot.vcf.gz"),
    message:
        "adding rough VAF to longshot variant calls for {wildcards.sample}"
    threads: 1
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/add_vaf.log.txt"),
    resources:
        cpus=1,
    shell:
        """
        sed -e '4i##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">' -e "s/SAMPLE/{wildcards.sample}/g" {input.vcf} | grep -v "DP\=0;" | \
        awk -v OFS="\t" -F"\t" '
        /^[^#]/{{ AC=$8; DP=$8;
        sub("DP=[0-9]*;", "", AC);
        sub("AC=[0-9]*,", "", AC);
        gsub(";.*", "", AC);
        sub(";.*", "", DP);
        sub("DP=", "", DP);
        $8 = $8"AF="AC/DP; }}1' | \
        bgzip -c > {output.vcf}
        """


# rule parse_medaka_vcf:
#     input:
#         vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.tmp.vcf"),
#     output:
#         vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf.gz"),
#     params:
#         parse_script = parse_script,
#         bam = os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
#         vcf = os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf"),
#         ref_index = config["reference"]+'.fai',
#     message:
#         "adding rough VAF to medaka variant calls for {wildcards.sample}"
#     threads: 1
#     log:
#         os.path.join(RESULT_DIR, "{sample}/logs/add_vaf.log.txt"),
#     resources:
#         cpus=1,
#     conda:
#         "../envs/medaka.yaml"
#     shell:
#         """
#         python3 {params.parse_script} -v {input.vcf} -b {params.bam} -r {params.ref_index} -o {params.vcf}
#         """


#rule parse_medaka_vcf:
#   input:
#       vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.tmp.vcf"),
#   output:
#       vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf.gz"),
#   message:
#       "adding rough VAF to medaka variant calls for {wildcards.sample}"
#   threads: 4
#   log:
#       os.path.join(RESULT_DIR, "{sample}/logs/add_vaf.log.txt"),
#   resources:
#       cpus=4,
#   shell:
#       """
#       sed -e '8i##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">' \
#           -e '9i##INFO=<ID=SAC,Number=1,Type=Integer,Description="Summed alt count">' \
#           -e '10i##INFO=<ID=NQ,Number=1,Type=Integer,Description="Normalised quality (QUAL/DP)">' \
#           -e "s/SAMPLE/{wildcards.sample}/g" {input.vcf} | \
#       grep -v "DPSP\=0;" | \
#       awk -v OFS="\t" -F"\t" '
#       /^[^#]/{{ AC=$8; DP=$8; DPSP=$8;
#				sub("AR=.*SR=", "SR=", AC);
#              sub("SR=[0-9]*,[0-9]*,", "", AC);
#				AC1=AC; AC2=AC;
#				sub(",[0-9]*","",AC1);
#				sub("[0-9]*,","",AC2);
#				AC=AC1+AC2;
#				sub("AR=.*DPSP=", "DPSP=", DPSP);
#				sub(";.*", "", DPSP);
#               sub("DPSP=", "", DPSP);
#				sub("AR=.*DP=", "DP=", DP);
#				sub(";.*", "", DP);
#               sub("DP=", "", DP);
#               $8 = $8";SAC="AC";AF="AC/DPSP";NQ="$6/DP; }}1' | \
#       bgzip -c > {output.vcf}
#       """

##### MEDAKA VARIANT CALLING END #####

##### CLAIR3 VARIANT CALLING START #####
rule clair3_variant:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        clair3_vcf=os.path.join(RESULT_DIR, "{sample}","clair3","merge_output.vcf.gz"),
    message:
        "clair3 variant calls for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/clair3.log.txt"),
    params:
        reference=config["reference"],
        refname=os.path.basename(config["reference"]).replace(".fasta",""),
        candidate_bed_path=os.path.join(RESULT_DIR, "{sample}","clair3","tmp/full_alignment_output/candidate_bed"),
        model=clair3_model,
        output=os.path.join(RESULT_DIR, "{sample}","clair3"),
    resources:
        cpus=4,
        #gpu=1,
    container:
        "docker://hkubal/clair3:latest"
    shell:
        """
        mkdir -p {params.candidate_bed_path}
        touch "{params.candidate_bed_path}/FULL_ALN_FILE_{params.refname}"
        /opt/bin/run_clair3.sh --bam_fn={input.bam} --sample_name={wildcards.sample} --ref_fn={params.reference} --threads={threads} --platform="ont" --model_path="/opt/models/{params.model}" --output={params.output}    --chunk_size=29903 --include_all_ctgs --no_phasing_for_fa --remove_intermediate_dir --enable_long_indel --haploid_sensitive 2&>{log}
        """


rule cleanup_clair3:
    input:
        clair3_vcf=os.path.join(RESULT_DIR, "{sample}","clair3","merge_output.vcf.gz"),
    output:
        clair3_vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.clair3.vcf.gz"),
    message:
        "reorganising clair3 variant calls for {wildcards.sample}"
    threads: 1
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/clair3.log.txt"),
    params:
        reference=config["reference"],
        model=clair3_model,
        output=os.path.join(RESULT_DIR, "{sample}","clair3")
    resources:
        cpus=1,
    shell:
        """
        bcftools view {input.clair3_vcf} | \
        sed -e '7i##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">' \
        -e '7i##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">' | \
        grep -v "DP\=0;" | \
        awk -v OFS="\t" -F"\t" '
        /^[^#]/{{ DP=$10; AF=$10;
        sub("[0-9]*:[0-9]*:", "", DP);
        sub(":.*", "", DP);
        sub(".*:","",AF)
        $8 = $8";DP="DP";AF="AF; }}1' | \
        bgzip -c > {output.clair3_vcf}
        """

##### CLAIR3 VARIANT CALLING END #####


##### LOFREQ VARIANT CALLING START #####
rule lofreq_dindel:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        bam=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.dindel.bam")),
    message:
        "lofreq dindel for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/lofreq.log.txt"),
    params:
        reference=config["reference"],
    resources:
        cpus=4,
        #gpu=1,
    conda:
        "../envs/lofreq.yaml"
    shell:
        """
        lofreq indelqual --dindel {input.bam} -f {params.reference} -o {output.bam} 2> /dev/null





        """

rule index_dindel:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.dindel.bam"),
    output:
        bam=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.dindel.bam.bai")),
        ckp=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.dindel_indexed.txt")),
    threads: 1
    resources:
        cpus=1,
        #gpu=1,
    shell:
        """
        samtools index -@{threads} {input.bam}
        touch {output.ckp}
        """

rule lofreq_call:
    input:
        ckp=os.path.join(RESULT_DIR, "{sample}/{sample}.dindel_indexed.txt"),
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.dindel.bam"),
    output:
        lofreq_vcf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.lofreq.vcf")),
    message:
        "lofreq variant calls for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/lofreq.log.txt"),
    params:
        reference=config["reference"],
    resources:
        cpus=4,
        #gpu=1,
    conda:
        "../envs/lofreq.yaml"
    shell:
        """
        lofreq call-parallel --no-baq --call-indels --pp-threads {threads} \
        -f {params.reference} -o {output.lofreq_vcf} {input.bam} 2> {log}
        """

rule cleanup_lofreq:
    input:
        lofreq_vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.lofreq.vcf"),
    output:
        lofreq_vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.lofreq.vcf.gz"),
    message:
        "reorganising lofreq variant calls for {wildcards.sample}"
    threads: 1
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/lofreq.log.txt"),
    params:
        reference=config["reference"],
    resources:
        cpus=1,
    shell:
        """
        sed -e '6i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
        -e "s|FILTER\tINFO|FILTER\tINFO\tFORMAT\t{wildcards.sample}|g" {input.lofreq_vcf} | \
        awk -F'\t' -v genotype=1 -v OFS="\t" '/^[^#]/{{ $9 = "GT"; $10 = genotype }}1' | \
        bgzip -c > {output.lofreq_vcf}
        bcftools index -f {output.lofreq_vcf}
        """

##### lofreq VARIANT CALLING END #####


rule filter_vcf:
    input:
        vcf_file=expand(os.path.join(RESULT_DIR, "{{sample}}/{{sample}}.{ext}"), ext=filter_extension),
    output:
        vcf_file=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.filtered.vcf.gz")),
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/{sample}.bcftools_filtering.log"),
    params:
        snv_freq=config["snv_min_freq"],
        snv_min_depth=config["min_depth"],
        snv_min_qual=snv_min_qual
    message:
        "setting conditional GT for {wildcards.sample}"
    shell:
        """
        bcftools index -f {input.vcf_file}
        bcftools +fill-tags {input.vcf_file} -Ou -- -t "TYPE" | \
        bcftools norm -Ou -a -m -  2> /dev/null | \
        bcftools view -f 'PASS,dn,dp,.' -i "INFO/AF >= {params.snv_freq} && INFO/DP >= {params.snv_min_depth} && QUAL >= {params.snv_min_depth}" -Oz -o {output.vcf_file}
        bcftools +setGT {output.vcf_file} -o {output.vcf_file} -- -t a -n 'c:1/1' 2>> {log}
        bcftools index {output.vcf_file}
        """


rule set_vcf_genotype:
    input:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/{sample}.filtered.vcf.gz"),
    output:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/{sample}.final.vcf.gz"),
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/{sample}.bcftools_setGT.log"),
    params:
        snv_freq=config["snv_min_freq"],
        con_freq=config["consensus_freq"],
        indel_freq=config["indel_freq"],
    message:
        "setting conditional GT for {wildcards.sample}"
    shell:
        """
        cp {input.vcf_file} {output.vcf_file}
        bcftools +setGT {output.vcf_file} -- -t q -i 'GT="1/1" && INFO/AF < {params.con_freq}' -n 'c:0/1' 2>> {log} | \
        bcftools +setGT -- -t q -i 'TYPE="indel" && INFO/AF < {params.indel_freq}' -n . 2>> {log} | \
        bcftools +setGT -o {output.vcf_file} -- -t q -i 'GT="1/1" && INFO/AF >= {params.con_freq}' -n 'c:1/1' 2>> {log}
        bcftools index -f {output.vcf_file}
        """


rule generate_consensus:
    input:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/{sample}.final.vcf.gz"),
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        consensus=os.path.join(RESULT_DIR, "{sample}/{sample}.consensus.fasta"),
        mask=temp(os.path.join(RESULT_DIR, "{sample}/mask.bed")),
        variants_bed=temp(os.path.join(RESULT_DIR, "{sample}/variants.bed")),
    message:
        "calling a consensus for {wildcards.sample}"
    threads: 1
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/{sample}.generate_consensus.log"),
    params:
        prefix="{sample}",
        reference=config["reference"],
        freq=config["consensus_freq"],
        min_depth=config["min_depth"],
    resources:
        cpus=1,
    shell:
        """
        bcftools query -f'%CHROM\t%POS0\t%END\n' {input.vcf_file} > {output.variants_bed}
        varCheck=$(file {output.variants_bed} | cut -f2 -d " ")
        if [ $varCheck == "empty" ]; then
            echo "BAM file was empty for {params.prefix}. Making empty consensus genome."
            emptySeq=$(printf %.1s N{{1..29903}})
            printf ">{params.prefix}\n$emptySeq\n" > {output.consensus}
            touch {output.mask}
            touch {output.variants_bed}
        else
            bedtools genomecov -bga -ibam {input.bam} | awk '$4 < {params.min_depth}' | \
            bedtools subtract -a - -b {output.variants_bed} > {output.mask}
            bcftools consensus -p {params.prefix} -f {params.reference} --mark-del '-' -m {output.mask} -H I -i 'INFO/DP >= {params.min_depth} & GT!="mis"' {input.vcf_file} 2> {log} | \
            sed "/^>/s/{params.prefix}.*/{params.prefix}/" > {output.consensus}
        fi
        """


rule amino_acid_consequences:
    input:
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.filtered.vcf.gz"),
    output:
        tsv=os.path.join(RESULT_DIR, "{sample}/{sample}.annotated.tsv"),
    message:
        "annotating variants for {wildcards.sample}"
    threads: 1
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/{sample}.aa_consequences.log"),
    params:
        reference=config["reference"],
        annotation=config["annotation"],
    resources:
        cpus=1,
    shell:
        """
        printf "reference\tsample_id\tgene\tnt_pos\tref\talt\teffect\taa_bcsq\taa_standard\tvaf\tdepth\n" > {output.tsv}
        bcftools csq -f {params.reference} -g {params.annotation} --force {input.vcf} 2>{log} | \
        bcftools query -f'[%CHROM\t%SAMPLE\t%POS\t%REF\t%ALT\t%DP\t%AF\t%TBCSQ\n]' 2>{log} | \
        awk -F'|' -v OFS="\t" '{{ print $1,$2,$6 }}' | \
        awk -v OFS="\t" -F"\t" '
        {{ POS=$10; REF=$10;
        gsub("[A-Z*].*", "", POS);
        gsub("[0-9]*", "", REF);
        ALT=REF;
        gsub(">.*", "", REF);
        gsub(".*>", "", ALT);
        NEW = sprintf("%s%s%s", REF, POS, ALT);
        $11 = NEW ; }} 1 ' | \
        awk -v OFS="\t" '{{ print $1,$2,$9,$3,$4,$5,$8,$10,$11,$7,$6 }}' >> {output.tsv}
        """

rule update_pangolin:
    output:
        update_info = os.path.join(RESULT_DIR, "pangolin_update_info.txt"),
    conda:
        "../envs/pangolin.yaml"
    shell:
        """
        echo "Attempting to update pangolin..." > {output.update_info}
        pangolin --update &>> {output.update_info} || echo "pangolin couldn't update. Investigate." >> {output.update_info}
        echo "Pangolin versions after attempted update..." >> {output.update_info}
        pangolin --all-versions &>> {output.update_info}
        """

rule pangolin:
    input:
        update_info = os.path.join(RESULT_DIR, "pangolin_update_info.txt"),
        fasta=os.path.join(RESULT_DIR, "{sample}/{sample}.consensus.fasta"),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/{sample}.lineage_report.csv"),
    conda:
        "../envs/pangolin.yaml"
    resources:
        cpus=1,
    threads: 4,
    shell:
        """
        pangolin --skip-scorpio --outfile {output.report} {input.fasta} &> /dev/null
        """


rule update_nextclade:
    output:
        update_info = os.path.join(RESULT_DIR, "nextclade_update_info.txt"),
    params:
        nextclade_dataset = config['nextclade_dataset']
    container:
        "docker://nextstrain/nextclade:2.3.1"
    shell:
        """
        echo "nextclade version:" > {output.update_info}
        nextclade --version >> {output.update_info} &>/dev/null
        echo "Updating SARS-CoV-2 dataset..." >> {output.update_info}
        nextclade dataset get --name sars-cov-2 -o {params.nextclade_dataset} &>>{output.update_info}
        """

rule nextclade:
    input:
        update_info = os.path.join(RESULT_DIR, "nextclade_update_info.txt"),
        fasta=os.path.join(RESULT_DIR, "{sample}/{sample}.consensus.fasta"),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/{sample}.nextclade_report.tsv"),
    params:
        nextclade_dataset = config['nextclade_dataset'],
        outdir = os.path.join(RESULT_DIR, "{sample}/nextclade"),
    container:
        "docker://nextstrain/nextclade:2.3.1"
    resources:
        cpus=1,
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/{sample}.nextclade.log"),
    threads: 4,
    shell:
        """
        nextclade run --in-order \
        --input-dataset={params.nextclade_dataset} \
        --output-basename={wildcards.sample} \
        --output-all={params.outdir} \
        --output-tsv={output.report} \
        --jobs 4 \
        {input.fasta} &> {log}
        """


#rule get_coverage_old:
#    input:
#        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
#    output:
#        bedtools_coverage=temp(
#            os.path.join(RESULT_DIR, "{sample}/{sample}.bedtools_coverage.tsv")
#        ),
#        samtools_coverage=temp(
#            os.path.join(RESULT_DIR, "{sample}/{sample}.samtools_coverage.tsv")
#        ),
#    resources:
#        cpus=1,
#    threads: 1
#    shell:
#        """
#        bedtools genomecov -ibam {input.bam} -d > {output.bedtools_coverage} 2>/dev/null
#        samtools coverage {input.bam} -o {output.samtools_coverage} 2>/dev/null
#        """

rule get_coverage:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        samtools_depth=temp(
            os.path.join(RESULT_DIR, "{sample}/{sample}.samtools_depth.tsv")
        ),
        samtools_coverage=temp(
            os.path.join(RESULT_DIR, "{sample}/{sample}.samtools_coverage.tsv")
        ),
    params:
        bed = config['coverage_bed'],
    resources:
        cpus=1,
    threads: 1
    shell:
        """
        samtools depth -a -J -@1 -b {params.bed} {input.bam} -o {output.samtools_depth}
        samtools coverage {input.bam} -o {output.samtools_coverage} 2>/dev/null
        """


rule sample_qc:
    input:
        samtools_depth=os.path.join(
            RESULT_DIR, "{sample}/{sample}.samtools_depth.tsv"
        ),
        samtools_coverage=os.path.join(
            RESULT_DIR, "{sample}/{sample}.samtools_coverage.tsv"
        ),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/{sample}.qc_results.csv"),
    params:
        fastq=os.path.join(config["reads_dir"], "{sample}.fastq"),
        sample="{sample}",
        min_depth=config["min_depth"],
    resources:
        cpus=1,
    threads: 1
    run:
        df = pd.read_csv(
            input.samtools_depth,
            sep="\t",
            header=None,
            names=["chrom", "pos", "depth"],
        )
        sam_df = pd.read_csv(input.samtools_coverage, sep="\t")
        num_mapped_reads = sam_df.loc[0, "numreads"]
        try:
            coverage_value = (df[df["depth"] >= params.min_depth].shape[0]) / (df.shape[0]) * 100
            mean_depth = df["depth"].mean()
        except:
            coverage_value = 0
            mean_depth = 0
        cmd = "echo $(cat {0}|wc -l)/4|bc".format(params.fastq)
        num_reads = int(os.popen(cmd).read().strip())

        out_df = pd.DataFrame(
            columns=[
                "id",
                "num_reads",
                "num_mapped_reads",
                coverage_colname,
                "mean_depth",
                "coverage_QC",
                "technology",
            ]
        )

        if coverage_value > 79:
            out_df.loc[0] = [
                params.sample,
                num_reads,
                num_mapped_reads,
                coverage_value,
                mean_depth,
                "PASS",
                "Oxford Nanopore MinION",
            ]
        else:
            out_df.loc[0] = [
                params.sample,
                num_reads,
                num_mapped_reads,
                coverage_value,
                mean_depth,
                "FAIL",
                "Oxford Nanopore MinION",
            ]

        out_df.to_csv(output.report, header=True, index=False)
