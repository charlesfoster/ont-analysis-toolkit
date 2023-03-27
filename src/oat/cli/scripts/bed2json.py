# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import sys
import re
import argparse
import json
import pandas as pd

## See https://stackoverflow.com/a/480227
def uniqify(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def main(sysargs=sys.argv[1:]):
    #set up argparse
    parser=argparse.ArgumentParser(description="bed2json.py: converts bedfile with amplicon primer locations into a primers.json file compatible with RAMPART", formatter_class=argparse.RawTextHelpFormatter, epilog="")
    parser.add_argument("bedfile", nargs="*", help="Path to bedfile with coordinates for where primers match the reference genome.")
    parser.add_argument('-l','--left_suffix', required=False, action='store', default = "_LEFT", help="Suffix for left primer in name column of BED file (default: '_LEFT').")
    parser.add_argument('-r','--right_suffix', required=False, action='store', default = "_RIGHT", help="Suffix for right primer in name column of BED file (default: '_RIGHT').")
    parser.add_argument('-o','--outfile', required=False, action='store', default = os.path.join(os.getcwd(),'primers.json'), help="Outfile name. Note that RAMPART require the file to (at least eventually) be called 'primers.json' and put into the correct directory.")

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    bedfile = ''.join(args.bedfile)
    if not os.path.isfile(bedfile):
        sys.exit("Error: bedfile does not appear to exist\n")

    df = pd.read_csv(bedfile,
                     sep = "\t",
                     header=None)
    
    amplicons = uniqify([re.sub('_LEFT|_RIGHT','',x) for x  in df.loc[:,3]])
    
    regions = []
    for amplicon in amplicons:
        start = min(list(df[df.loc[:,3].str.contains(amplicon, case=False, na=False)].loc[:,1]))
        end = max(list(df[df.loc[:,3].str.contains(amplicon, case=False, na=False)].loc[:,2]))
        regions.append([start,end])
    
    out_dict = {'name':'Primer scheme',
                'amplicons':regions}
    
    # write to file
    fname = os.path.join(args.outfile)
    with open(fname, "w") as fp:
        json.dump(out_dict, fp, indent=2)
    
    print("\nConverted {0} to {1}".format(bedfile, args.outfile))

if __name__ == '__main__':
    main()
    