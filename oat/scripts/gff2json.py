#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 11:45:20 2022

@author: cfos
"""

from Bio import SeqIO
import json
import os
import sys
import argparse
import pandas as pd
import re

def main(sysargs=sys.argv[1:]):
    #set up argparse
    parser=argparse.ArgumentParser(description="gff2json.py: converts gff3-format file + fasta-format reference genome into a genome.json file compatible with RAMPART", formatter_class=argparse.RawTextHelpFormatter, epilog="")
    parser.add_argument('-f',"--fasta", required=True, help="Path to fasta file of reference genome.")
    parser.add_argument('-g',"--gff3", required=True, help="Path to gff3 format annotation file for reference genome.")
    parser.add_argument('-o','--outfile', required=False, action='store', default = os.path.join(os.path.abspath(os.path.dirname(__file__)),'genome.json'), help="Outfile name. Note that RAMPART require the file to (at least eventually) be called 'genome.json' and put into the correct diectory.")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    gff3_file = ''.join(args.gff3)   
    if not os.path.isfile(gff3_file):
        sys.exit("Error: gff3 file does not appear to exist\n")

    fasta_file = ''.join(args.fasta)   
    if not os.path.isfile(fasta_file):
        sys.exit("Error: fasta file does not appear to exist\n")
    

    # read in and wrangle gff3 file
    df = pd.read_csv(gff3_file,
                     comment = "#",
                     sep = "\t",
                     header=None,
                     names=['seqid','source','type','start','end','score','strand','phase','attributes'])
    
    length = int(df.iloc[0]['end'])
    label = str(df.iloc[0]['seqid'])
    df = df[df['type'] == 'gene']
    
    genes_list = [re.search('Name\=(.+?);', x).group(1) for x in list(df['attributes'])]
    df['attributes'] = genes_list
    dupes = list(set(sorted(genes_list)[::2]) & set(sorted(genes_list)[1::2]))
    
    # get genes part of json
    genes_dict = {}
    for gene in genes_list:
        if gene not in dupes:
            strand = 1 if str(df[df['attributes']==gene]['strand'].iat[0]) == '+' else -1
            gene_dict = {'start':int(df[df['attributes']==gene]['start']),
                         'end':int(df[df['attributes']==gene]['end']),
                         'strand':strand}
            genes_dict[gene] = gene_dict
        elif gene in dupes:
            sub_df = df[df['attributes']==gene]
            for i in range(len(sub_df)):
                name = gene+f'_anno{i+1}'
                strand = '1' if str(sub_df['strand'].iat[0]) == '+' else '-1'
                gene_dict = {'start':int(sub_df.iloc[i]['start']),
                             'end':int(sub_df.iloc[i]['end']),
                             'strand':strand}
                genes_dict[name] = gene_dict
    
    # get reference part of json
    fasta = next(SeqIO.parse(fasta_file, "fasta"))
    
    reference = {'label':label,
                 'accession':df.iat[0,0],
                 'sequence':str(fasta.seq)}
    
    # put it all together
    final = {'label':label,
             'length':length,
             'genes':genes_dict,
             'reference':reference}
    
    # write to file
    fname = os.path.join(args.outfile)
    with open(fname, "w") as fp:
        json.dump(final, fp, indent=2)
    
    print("\nConverted {0} to {1}".format(gff3_file, args.outfile))

if __name__ == '__main__':
    main()
    
