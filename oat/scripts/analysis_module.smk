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
# set up conda path for pangolin
cmd = 'eval "$(conda shell.bash hook)" && echo $(conda info --base)'
conda_path = os.popen(cmd).read()
conda_root = os.path.split(conda_path)[0]

# set up env variable & threads for tensorflow / medaka
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
medaka_con_threads = int(config["threads"] / 2)

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
        if os.path.getsize(file) == 0 or file.endswith(
            (
                ".all.vcf.gz.csi",
                "draft.vcf.gz.tbi",
                ".filtered.vcf.gz.csi",
                ".mapped.bam.csi",
            )
        ):
            os.remove(file)


################
# rules
################

potential_pangolin = []
if SARS_ANALYSIS:
    potential_pangolin.append(
        expand(
            os.path.join(RESULT_DIR, "{sample}/{sample}.lineage_report.csv"),
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
        fa_files = glob.glob(RESULT_DIR + "/**/*.consensus.fasta", recursive=True)

        multifasta = os.path.join(
            RESULT_DIR, config["run_name"] + ".consensus_genomes.fa"
        )
        if os.path.exists(multifasta):
            os.remove(multifasta)

        for fa in fa_files:
            os.system("cat {0} >> {1}".format(fa, multifasta))
        qc_files = glob.glob(RESULT_DIR + "/**/*qc_results.csv", recursive=True)
        combined_qc = pd.concat([pd.read_csv(f) for f in qc_files]).set_index(["id"])
        run_metadata = pd.read_csv(params.run_metadata)
        outdata = run_metadata.join(combined_qc, on=["id"])
        outdata = outdata.set_index("id", drop=False)
        outdata.index.name = None

        # add analysis metadata
        outdata["analysis_date"] = date.today().strftime("%Y-%m-%d")

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

        # define order of columns
        if params.sars_analysis:
            compulsory_col_order = [
                "analysis_date",
                "run_name",
                "id",
                "barcode",
                "lineage",
                "scorpio_call",
                "ref_cov_20",
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
                "ref_cov_20",
                "mean_depth",
                "num_reads",
                "num_mapped_reads",
                "percent_total_reads",
                "reads_qc",
            ]

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


rule map_reads:
    input:
        fastq=os.path.join(config["reads_dir"], "{sample}.fastq"),
    output:
        bam=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.mapped.bam")),
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
        samtools ampliconclip --both-ends --soft-clip --filter-len 30 -@ 20 -u -b {params.bed} {input.bam} 2>{log} | \
        samtools sort -u -@ 20 -n 2>>{log} | \
        samtools fixmate -u -@ 20 - - 2>>{log} | \
        samtools calmd -u -@ 20 - {params.reference} 2>>{log} | \
        samtools sort -u -@ 20 2>>{log} | \
        samtools view --write-index -u -@ 20 -F 4 -o {output.bam}
        """


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
        model="r941_min_high_g360",
    resources:
        cpus=medaka_con_threads,
    conda:
        "../envs/medaka.yaml"
    shell:
        """
        medaka consensus --model {params.model} --threads {threads} --chunk_len 800 --chunk_ovlp 400 {input.bam} {output.hdf} 2>{log}
        """


rule medaka_variant:
    input:
        hdf=os.path.join(RESULT_DIR, "{sample}/{sample}.hdf"),
    output:
        medaka_vcf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.medaka.vcf")),
        vcf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.draft.vcf.gz")),
    message:
        "initial medaka variant calls for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/medaka_variant.log.txt"),
    params:
        reference=config["reference"],
        model="r941_min_high_g360",
    resources:
        cpus=4,
    conda:
        "../envs/medaka.yaml"
    shell:
        """
        medaka variant {params.reference} {input.hdf} {output.medaka_vcf} 2>{log}
        bgzip -c {output.medaka_vcf} > {output.vcf} 2>>{log}
        tabix -p vcf {output.vcf} 2>>{log}
        """


rule longshot:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.draft.vcf.gz"),
    output:
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
        # longshot -P 0 -F -A --no_haps --bam {input.bam} --ref {params.reference} --out {output.vcf} --potential_variants {input.vcf}    2>{log}
         longshot -P 0 -F --no_haps --bam {input.bam} --ref {params.reference} --out {output.vcf} --potential_variants {input.vcf}    2>{log}
        """


rule add_rough_VAF:
    input:
        vcf=os.path.join(RESULT_DIR, "{sample}/{sample}.all.vcf"),
    output:
        vcf=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.all.vcf.gz")),
    message:
        "adding rough VAF to longshot variant calls for {wildcards.sample}"
    threads: 4
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/add_vaf.log.txt"),
    resources:
        cpus=4,
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


rule filter_vcf:
    input:
        vcf_file=os.path.join(RESULT_DIR, "{sample}/{sample}.all.vcf.gz"),
    output:
        vcf_file=temp(os.path.join(RESULT_DIR, "{sample}/{sample}.filtered.vcf.gz")),
    log:
        os.path.join(RESULT_DIR, "{sample}/logs/{sample}.bcftools_filtering.log"),
    params:
        snv_freq=config["snv_min"],
        con_freq=config["consensus_freq"],
        snv_min=20,
    message:
        "setting conditional GT for {wildcards.sample}"
    shell:
        """
        bcftools index -f {input.vcf_file}
        bcftools +fill-tags {input.vcf_file} -Ou -- -t "TYPE" | \
        bcftools norm -Ou -a -m -  2> /dev/null | \
        bcftools view -f 'PASS,dn,dp,.' -i "INFO/AF >= {params.snv_freq} && INFO/DP >= {params.snv_min}" -Oz -o {output.vcf_file}
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
        snv_freq=config["snv_min"],
        con_freq=config["consensus_freq"],
    message:
        "setting conditional GT for {wildcards.sample}"
    shell:
        """
        cp {input.vcf_file} {output.vcf_file} 
        bcftools +setGT {output.vcf_file} -- -t q -i 'GT="1/1" && INFO/AF < {params.con_freq}' -n 'c:0/1' 2>> {log} | \
        bcftools +setGT -- -t q -i 'TYPE="indel" && INFO/AF < {params.con_freq}' -n . 2>> {log} | \
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
            bedtools genomecov -bga -ibam {input.bam} | awk '$4 < 20' | \
            bedtools subtract -a - -b {output.variants_bed} > {output.mask}
            bcftools consensus -p {params.prefix} -f {params.reference} --mark-del '-' -m {output.mask} -H I -i 'INFO/DP >= 20 & GT!="mis"' {input.vcf_file} 2> {log} | \
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


rule pangolin:
    input:
        fasta=os.path.join(RESULT_DIR, "{sample}/{sample}.consensus.fasta"),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/{sample}.lineage_report.csv"),
    params:
        pangolin_path=os.path.join(conda_root, "pangolin"),
    shell:
        """
        set +eu
        #eval "$(conda shell.bash hook)" && conda activate {params.pangolin_path} && pangolin --outfile {output.report} {input.fasta} &> /dev/null
        eval "$(conda shell.bash hook)" && conda activate pangolin && pangolin --outfile {output.report} {input.fasta} &> /dev/null
        set -eu
        """


rule get_coverage:
    input:
        bam=os.path.join(RESULT_DIR, "{sample}/{sample}.trimmed.bam"),
    output:
        bedtools_coverage=temp(
            os.path.join(RESULT_DIR, "{sample}/{sample}.bedtools_coverage.tsv")
        ),
        samtools_coverage=temp(
            os.path.join(RESULT_DIR, "{sample}/{sample}.samtools_coverage.tsv")
        ),
    resources:
        cpus=1,
    threads: 1
    shell:
        """
        bedtools genomecov -ibam {input.bam} -d > {output.bedtools_coverage} 2>/dev/null
        samtools coverage {input.bam} -o {output.samtools_coverage} 2>/dev/null
        """


if SARS_ANALYSIS:
    sample_qc_lineages = os.path.join(
        RESULT_DIR, "{sample}/{sample}.lineage_report.csv"
    )
else:
    sample_qc_lineages = ([None],)


rule sample_qc:
    input:
        lineages=sample_qc_lineages,
        bedtools_coverage=os.path.join(
            RESULT_DIR, "{sample}/{sample}.bedtools_coverage.tsv"
        ),
        samtools_coverage=os.path.join(
            RESULT_DIR, "{sample}/{sample}.samtools_coverage.tsv"
        ),
    output:
        report=os.path.join(RESULT_DIR, "{sample}/{sample}.qc_results.csv"),
    params:
        sars_analysis=SARS_ANALYSIS,
        fastq=os.path.join(config["reads_dir"], "{sample}.fastq"),
        sample="{sample}",
    resources:
        cpus=1,
    threads: 1
    run:
        df = pd.read_csv(
            input.bedtools_coverage,
            sep="\t",
            header=None,
            names=["chrom", "pos", "depth"],
        )
        sam_df = pd.read_csv(input.samtools_coverage, sep="\t")
        num_mapped_reads = sam_df.loc[0, "numreads"]
        ref_cov_20 = (df[df["depth"] >= 20].shape[0]) / (df.shape[0]) * 100
        mean_depth = df["depth"].mean()
        cmd = "echo $(cat {0}|wc -l)/4|bc".format(params.fastq)
        num_reads = int(os.popen(cmd).read().strip())

        out_df = pd.DataFrame(
            columns=[
                "id",
                "num_reads",
                "num_mapped_reads",
                "ref_cov_20",
                "mean_depth",
                "coverage_QC",
            ]
        )

        if ref_cov_20 > 79:
            out_df.loc[0] = [
                params.sample,
                num_reads,
                num_mapped_reads,
                ref_cov_20,
                mean_depth,
                "PASS",
            ]
        else:
            out_df.loc[0] = [
                params.sample,
                num_reads,
                num_mapped_reads,
                ref_cov_20,
                mean_depth,
                "FAIL",
            ]

        if params.sars_analysis:

            lineages = pd.read_csv(input.lineages).drop(
                ["conflict", "note", "scorpio_conflict"], axis=1
            )
            lineages.columns = [
                "id",
                "lineage",
                "pangolin_ambiguity_score",
                "scorpio_call",
                "scorpio_support",
                "lineage_designation_version",
                "pangolin_version",
                "pangoLEARN_version",
                "pango_version",
                "pangolin_status",
            ]
            lineages.loc[0, "id"] = wildcards.sample

            out_df = out_df.join(lineages.set_index(["id"]), on=["id"])

        out_df.to_csv(output.report, header=True, index=False)
