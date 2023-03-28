#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 09:48:54 2022

@author: cfos
"""

from cyvcf2 import VCF, Writer
from medaka.features import pileup_counts
from medaka.common import Region
import numpy as np
import statistics
import argparse
import sys 
import os

#==============================================================#
# %% DEFINE ACCESSORY FUNCTIONS                                #
#==============================================================#

def base_lookup(alt_base, counts_matrix, position):
    #print(alt_base)
    #print(position)
    array = counts_matrix[position][0]
    
    counts = {'A':array[0] + array[4],
              'C':array[1] + array[5],
              'G':array[2] + array[6],
              'T':array[3] + array[7],
              'del':array[8] + array[9]}
    
    #print(counts)
        
    #print(f"{base}: {count}")
    result = {'VC':counts[alt_base], 'RDP':np.sum(array),'AF':counts[alt_base]/np.sum(array)}
    return(result)
        
def query_variant_snp(vcf_pos, alt_base, medaka_positions, medaka_counts):
    position = np.where((medaka_positions['major'] == vcf_pos-1) & (medaka_positions['minor'] == 0))
    result = base_lookup(alt_base, medaka_counts, position)
    return(result)

def query_variant_deletion(vcf_pos_start, vcf_pos_end, alt_base, medaka_positions, medaka_counts):
    RDP_list = []
    VC_list = []
    AF_list = []
    
    for pos in range(vcf_pos_start, vcf_pos_end):
        position = np.where((medaka_positions['major'] == pos-1) & (medaka_positions['minor'] == 0))
        result = base_lookup(alt_base, medaka_counts, position)
        RDP_list.append(int(result['RDP']))
        VC_list.append(int(result['VC']))
        AF_list.append(float(result['AF']))
    
    ave_result = {'RDP':round(statistics.mean(RDP_list)+0.1),
                  'VC':round(statistics.mean(VC_list)+0.1),
                  'AF':statistics.fmean(AF_list)}
    
    return(ave_result)

def query_variant_insertion(vcf_pos_start, vcf_pos_end, alt_string, medaka_positions, medaka_counts):
    RDP_list = []
    VC_list = []
    AF_list = []

    ctr = 1
    for pos in range(vcf_pos_start, vcf_pos_end-1):
        alt_base = alt_string[ctr]
        position = np.where((medaka_positions['major'] == vcf_pos_start) & (medaka_positions['minor'] == ctr))
        result = base_lookup(alt_base, medaka_counts, position)
        RDP_list.append(int(result['RDP']))
        VC_list.append(int(result['VC']))
        AF_list.append(float(result['AF']))
        ctr+=1
    
    ave_result = {'RDP':round(statistics.mean(RDP_list)+0.1),
                  'VC':round(statistics.mean(VC_list)+0.1),
                  'AF':statistics.fmean(AF_list)}
    
    return(ave_result)


#==============================================================#
# %% MAIN                                                      #
#==============================================================#

def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description="Use the medaka API to get base counts & frequencies for variants",
        usage="""python3 parse_medaka_variants.py [options] """,
        epilog="""
        This approach is very rough and should not be treated as rock-solid evidence. 
        The purpose is to give a rough idea of the underlying number of reads that
        support a variant, as well as the allele frequency (count/depth). Having
        these data can be useful in cases where a given sample has unexpected results, 
        e.g. unassigned lineage, and can reveal sites that have (e.g.) a roughly
        50/50 split of reads supporting different bases, which can _sometimes_ 
        indicate a QC issue or a mixed infection.
        
        The read depths and ALT counts for simple SNPs is calculated fairly
        easily, but these metrics are more difficult to calculate for indels using
        the approach in this script. As a compromise, the final metrics for indels 
        are averaged across the span of the variant.
        
        However, there are problems associated with purely taking the raw base counts as 
        evidence for variants, especially with Nanopore data. For example, 
        "counts of bases and evidence for variants are not synonymous".
        
        See the following posts:
        - https://community.nanoporetech.com/posts/calculating-variant-allele
        - https://community.nanoporetech.com/posts/what-is-the-best-approach
        - https://labs.epi2me.io/notebooks/Introduction_to_how_ONT's_medaka_works.html
        
        Again, treat results from this script as a rough guide for further investigation.
        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-b",
        "--bam",
        action="store",
        required=True,
        help="bam file for sample",
    )
    parser.add_argument(
        "-r",
        "--reference_index",
        action="store",
        required=True,
        help="Path to fai index of reference genome used for variant calling",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        action="store",
        required=False,
        default='edited.vcf',
        help="Output vcf name. Default: /path/to/input_vcf_directory/edited.vcf",
    )
    parser.add_argument(
        "-v",
        "--vcf",
        action="store",
        required=True,
        help="vcf file for sample",
    )

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    args.outfile = os.path.join(os.path.dirname(args.vcf),args.outfile)
    
    with open(args.reference_index, 'r') as f:
        ref_info = f.readlines()[0].split()

    region = Region.from_string(f'{ref_info[0]}:0-{ref_info[1]}')
    bam_file = args.bam
    pileup_data = pileup_counts(region, bam_file)
    pileup_data = pileup_data[0]  # implementation detail that need not trouble us
    medaka_counts, medaka_positions = pileup_data
    
    ### STEP 2: parse VCF and add in values
    VCF_IN = args.vcf
    VCF_OUT = args.outfile
    
    vcf = VCF(VCF_IN)
    vcf.add_info_to_header({'ID': 'RDP', 'Description': 'Raw count of reads spanning variant position calculation. For indels = mean across variant span.',
        'Type':'Integer', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'VC', 'Description': 'Count of reads supporting ALT. For indels = mean across variant span.',
        'Type':'Integer', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'AF', 'Description': 'Allele frequency of variant. For indels = mean across variant span.',
        'Type':'Float', 'Number': '1'})
    
    sample_name = os.path.basename(VCF_IN)
        
    w = Writer(VCF_OUT, vcf)
    
    for v in vcf:
        ref_base = v.REF
        alt_base = v.ALT[0]
        vcf_pos = v.POS
        if len(ref_base) == 1 and len(alt_base) == 1:
            try:
                results = query_variant_snp(vcf_pos, alt_base, medaka_positions, medaka_counts)
                v.INFO["RDP"] = int(results['RDP'])
                v.INFO["VC"] = int(results['VC'])
                v.INFO["AF"] = float(results['AF'])
                w.write_record(v)
                continue
            except:
                print("Skipping: {}".format(sample_name+":\t"+ref_base+str(vcf_pos)+alt_base))
                continue
        
        if len(ref_base) > 1 and len(alt_base) == 1:
            vcf_pos_start = vcf_pos+1
            if len(ref_base) == 2:
                vcf_pos_end = vcf_pos + len(ref_base) +1
            else:
                vcf_pos_end = vcf_pos + len(ref_base) 
            try:            
                results = query_variant_deletion(vcf_pos_start, vcf_pos_end, "del", medaka_positions, medaka_counts)
                v.INFO["RDP"] = int(results['RDP'])
                v.INFO["VC"] = int(results['VC'])
                v.INFO["AF"] = float(results['AF'])
                w.write_record(v)
                continue
            except:
                print("Skipping: {}".format(sample_name+":\t"+ref_base+str(vcf_pos)+alt_base))
                continue
    
        if len(ref_base) == 1 and len(alt_base) > 1:
            vcf_pos_start = vcf_pos-1
            if len(alt_base) == 2:
                vcf_pos_end = vcf_pos + len(alt_base) 
            else:
                vcf_pos_end = vcf_pos + len(alt_base) -1
            
            alt_string = alt_base[0:]
            try:
                results = query_variant_insertion(vcf_pos_start, vcf_pos_end, alt_string, medaka_positions, medaka_counts)
                v.INFO["RDP"] = int(results['RDP'])
                v.INFO["VC"] = int(results['VC'])
                v.INFO["AF"] = float(results['AF'])
                w.write_record(v)
                continue
            except:
                print("Skipping: {}".format(sample_name+":\t"+ref_base+str(vcf_pos)+alt_base))
                continue
    
    w.close(); vcf.close()
    os.system(f"bgzip -c {VCF_OUT} > {VCF_OUT+'.gz'}")
    os.system(f"tabix -p vcf {VCF_OUT+'.gz'}")

# %%                                                      #
if __name__ == '__main__':
    main()
    
