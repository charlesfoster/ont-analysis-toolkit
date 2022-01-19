#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 11:45:20 2022

@author: cfos
"""

from Bio import SeqIO
from Bio import SeqFeature
import json
import os
import sys
import argparse

def main(sysargs=sys.argv[1:]):
    #set up argparse
    parser=argparse.ArgumentParser(description="genbank2json.py: converts genbank format file into a genome.json file compatible with RAMPART", formatter_class=argparse.RawTextHelpFormatter, epilog="")
    parser.add_argument("genbank_file", nargs="*", help="Path to genbank format file for reference genome.")
    parser.add_argument('-o','--outfile', required=False, action='store', default = os.path.join(os.path.abspath(os.path.dirname(__file__)),'genome.json'), help="Outfile name. Note that RAMPART require the file to (at least eventually) be called 'genome.json' and put into the correct diectory.")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    
    gb_file = ''.join(args.genbank_file)
    
    if not os.path.isfile(gb_file):
        sys.exit("Error: genbank file does not appear to exist\n")
    
    record = next(SeqIO.parse(gb_file, "genbank"))
    
    # get genes part of json
    genes_dict = {}
    for f in record.features:
        if f.type == "gene":
            if isinstance(f.location, SeqFeature.CompoundLocation):
                for i in range(len(f.location.parts)):
                    name = ''.join(f.qualifiers['gene'])+'_part{}'.format(str(i+1))
                    gene_dict = {'start':int(f.location.parts[i].start)+1,
                                 'end':int(f.location.parts[i].end)+1,
                                 'strand':int(f.location.parts[i].strand)}
                    genes_dict[name] = gene_dict
    
            else:
                name = ''.join(f.qualifiers['gene'])
                gene_dict = {'start':int(f.location.start)+1,
                             'end':int(f.location.end)+1,
                             'strand':int(f.location.strand)}
                genes_dict[name] = gene_dict
    
    # get reference part of json
    reference = {'label':record.description,
                 'accession':record.id,
                 'sequence':str(record.seq)}
    
    # put it all together
    final = {'label':record.description,
             'length':len(record.seq),
             'genes':genes_dict,
             'reference':reference}
    
    # write to file
    fname = os.path.join(args.outfile)
    with open(fname, "w") as fp:
        json.dump(final, fp, indent=2)
    
    print("\nConverted {0} to {1}".format(gb_file, args.outfile))

if __name__ == '__main__':
    main()
    
