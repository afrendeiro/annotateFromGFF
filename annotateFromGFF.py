#!/usr/bin/env python
"""
annotateFromGFF.py

Creates functional annotation of complete genome based on features on GFF file.
Reports coding-sequence (CDS) and untranslated-regions (UTR) and extracts 
features such as introns, intergenic space and promoters. It also distinguishes
between 5' and 3' UTRs.
Promoters and intergenic space are annotated dynamically based on a specified
promoter size but rezisable to fit genome boundaries and the occurance of operons.

Andre F. Rendeiro (andre.rendeiro@mail.com) 2014

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.
"""

import sys, csv, re
from argparse import ArgumentParser


def main():
    # option parser    
    parser = ArgumentParser(description = 'functional annotation of complete genome based on features on GFF file.',
        usage = 'python annotateFromGFF.py [OPTIONS] file.gff chrmSizes.tsv > annotation.bed')
    # positional arguments
    parser.add_argument('gff',
        help = 'GFF file with annotation.')
    parser.add_argument('chrmFile',
        help = 'Tab-delimited file with sizes of each chromossome (chr:size).')
    parser.add_argument('-o', '--outfile', default = "", dest = 'outFile',
        help = 'Specify the name of the output file. If not specified, will output to stdout.')
    parser.add_argument('-p', '--promoterSize', type = int, dest = 'promSize', default = 300,
        help = 'Average size of promoter elements. Dynamically resizable.')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose',
        help = 'Print debugging information', default = False)
    # parse
    args = parser.parse_args()
    
    v = args.verbose

    # Start working

    chrmSizes = getChrSizes(args.chrmFile)
    output = parseGFF(args.gff, chrmSizes, args.promSize)
    writeOutput(output, args.outFile)

def getChrSizes(chrmFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    try:
        with open(chrmFile, 'r') as f:
            chrmSizes = {}
            for line in enumerate(f):
                row = line[1].strip().split('\t')
                chrmSizes[str(row[0])] = int(row[1])
        return chrmSizes
    except IOError:
        print("%s file not openable or doesn't exist" % chrmFile)
        sys.exit(0)

def parseGFF(infile, chrmSizes, promSize):
    """
    Parses GFF file and annotates genomic features.
    Promoters 
    Outputs list of lists with annotation for each feature.
    """
    try:
        with open(infile, 'r') as f:

            prev_line = ""
            prev_chrm = ""
            prev_end = 0
            cur_line = ""
            trigger = False

            output = []

            reader = csv.reader(f, delimiter='\t')

            for row in reader:

                chrm = str(row[0])
                start = int(row[3])
                end = int(row[4])
                cur_line = str(row[2])

                # start parsing gff
                if cur_line == "gene":
                    strand = str(row[6])
                    gene = re.search('GSOIDG\w+', str(row[8])).group(0)

                    if trigger:
                        istart = prev_end + 1
                        # check it's within chr boundaries
                        if prev_end + 1 + promSize <= chrmSizes[prev_chrm]:
                            iend = prev_end + 1 + promSize
                        else:
                            iend = chrmSizes[prev_chrm]
                        # Add TSS
                        cur_loc = "TSS"
                        output.append([prev_chrm, iend - 1, iend, cur_loc, gene])
                        # Add promoter from previous gene                        
                        cur_loc = "Promoter"
                        output.append([prev_chrm, istart, iend, cur_loc, gene])

                    # Add intergenic space
                    cur_loc = "Intergenic"

                    # check it's same chrm as before, if not fill rest of chromossome with
                    if (prev_chrm != chrm) and (prev_chrm != ""):
                        # fill end of last chromossomem with intergenic space 
                        if trigger:
                            istart = prev_end + 1 + 1 + promSize # 
                        else:
                            istart = prev_end + 1
                        iend = chrmSizes[prev_chrm]
                        output.append([prev_chrm, istart, iend, cur_loc])

                        # reset chrom
                        trigger = False
                        prev_end = 0

                    # check first if theres space for intergenic space
                    # get intergenic start
                    if trigger:
                        # check it's within chr boundaries
                        if prev_end + 1 + promSize >= 1:
                            istart = prev_end + 1 + promSize
                        else:
                            istart = 1
                    else:
                        # check it's within chr boundaries
                        if prev_end + 1 >= 1:
                            istart = prev_end + 1
                        else:
                            istart = 1
                    # get intergenic end
                    if strand == "+":
                        # check it's within chr boundaries
                        if start - 1 - promSize:
                            iend = start - 1 - promSize
                        else:
                            iend = 1
                    elif strand == "-":
                        # check it's within chr boundaries
                        if start - 1:
                            iend = start - 1
                        else:
                            iend = 1

                    # if it's enough add it
                    if iend - istart >= 1:
                        output.append([chrm, istart, iend, cur_loc])

                    trigger = False
                    #output.append([chrm, start, end, "Gene", strand, id])
                    
                    if strand == "+":
                        # Add promoter from current gene if positive
                        cur_loc = "Promoter"
                        # check there's space for promoter in chr
                        if start - 1 - promSize >= 1:
                            istart = start - 1 - promSize # if there's space
                        else:
                            istart = 1 # if there's no space reduce promoter size
                        if start - 1 >= 1:
                            iend = start - 1
                        else:
                            iend = 1
                        output.append([chrm, istart, iend, cur_loc, gene])
                        # Add TSS
                        cur_loc = "TSS"
                        output.append([chrm, istart, istart + 1, cur_loc, gene])

                    elif strand == "-":
                        trigger = True

                    start = istart
                    end = iend

                elif cur_line == "mRNA":
                    pass

                elif cur_line == "CDS":
                    if (prev_loc == "CDS") or ("UTR" in prev_loc):
                        # annotate intron
                        cur_loc = "Intron"
                        istart = prev_end + 1
                        iend = start - 1
                        output.append([chrm, istart, iend, cur_loc, gene])

                    # annotate CDS
                    cur_loc = "CDS"
                    output.append([chrm, start, end, cur_loc, gene])

                elif cur_line == "UTR":
                    # annotate intron
                    if (prev_loc == "CDS") or ("UTR" in prev_loc):
                        cur_loc = "Intron"
                        istart = prev_end + 1
                        iend = start - 1
                        output.append([chrm, istart, iend, cur_loc, gene])

                    # annotate 5' UTR
                    if (strand == "+" and prev_loc == "Intergenic") or (strand == "-" and (prev_loc == "CDS")):
                        cur_loc = "5'UTR"
                        output.append([chrm, start, end, cur_loc, gene])

                    # annotate 3'
                    elif (strand == "-" and prev_loc == "Intergenic") or (strand == "+" and (prev_loc == "CDS")):
                        cur_loc = "3'UTR"
                        output.append([chrm, start, end, cur_loc, gene])

                prev_line = cur_line
                prev_chrm = chrm
                prev_end = end
                prev_loc = cur_loc

        # Annotate chrsm without genes
        allChrms = set(chrmSizes.keys())
        curChrs = set([c[0] for c in output])

        for c in allChrms.difference():
            output.append([c, 1, chrmSizes[c], "Intergenic"])

        return output#.sort()

    except IOError:
        print(" %s file not openable or doesn't exist" % infile)
        sys.exit(0)

def writeOutput(output, outFile):
    """
    Writes output (list) tab-delimited
    to outFile or to stout if outFile is not specified.
    """
    try:
        if outFile != "":
            # write to file
            with open(outFile, 'wb') as f:
                wr = csv.writer(f, delimiter='\t')
                for line in output:
                    wr.writerow(line)
        else:
            wr = csv.writer(sys.stdout, delimiter='\t')
            for line in output:
                wr.writerow(line)
    except IOError:
        print(" %s file not writable" % outFile)
        sys.exit(0)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)