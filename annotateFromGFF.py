#!/usr/bin/env python
"""
annotateFromGFF.py

Creates functional annotation of complete genome based on features on GFF file.
Reports coding-sequence (CDS) and untranslated-regions (UTR) and extracts 
features such as introns, intergenic space, TSSs and promoters. It also distinguishes
between 5' and 3' UTRs.
Promoters and intergenic space are annotated dynamically based on a specified
promoter size but rezisable to fit genome boundaries and small intergenic space.

Andre F. Rendeiro (andre.rendeiro@mail.com) 2014

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.
"""

import sys, logging, csv, re
from argparse import ArgumentParser


def main():
    # argparser    
    parser = ArgumentParser(description = 'functional annotation of complete genome based on features on GFF file.',
        usage = 'python annotateFromGFF.py [OPTIONS] file.gff chrmSizes.tsv > annotation.bed')
    # positional arguments
    parser.add_argument('gff',
        help = 'GFF file with annotation.')
    parser.add_argument('chrmFile',
        help = 'Tab-delimited file with sizes of each chromossome (chr:size).')
    # optional arguments
    parser.add_argument('-o', '--outfile', default = "", dest = 'outFile',
        help = 'Specify the name of the output file. If not specified, will output to stdout.')
    parser.add_argument('-p', '--promoterSize', type = int, dest = 'promSize', default = 300,
        help = 'Average size of promoter elements. Dynamically resizable.')
    parser.add_argument('-op', '--operons', action = 'store_false', dest = 'operons',
        help = "Consider operons.", default = True) # default is "true for storing false" meaning off
    parser.add_argument('--operonDistance', type = int, dest = 'operonDist',
        help = "Distance to consider genes as belonging to same operon.", default = 60)
    parser.add_argument('-l', '--logfile', default = "log.txt", dest = 'logFile',
        help = 'Specify the name of the log file.')
    parser.add_argument('-s', '--silent', action = 'store_true', dest = 'silent',
        help = "Silent behaviour. Don't make log file.", default = False)
    # parse
    args = parser.parse_args()

    # logging
    global logger
    logger = logging.getLogger(sys.argv[0])
    logger.setLevel(logging.INFO)

    # create a file handler
    handler = logging.FileHandler(args.logFile)
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    global s
    s = args.silent

    # Start working
    chrmSizes = getChrSizes(args.chrmFile)
    output = parseGFF(args.gff, chrmSizes, args.promSize, args.operons, args.operonDist)
    writeOutput(output, args.outFile)
    if not s:
        logger.info("Finished run successfully.")
    sys.exit(1)

def getChrSizes(chrmFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    if not s:
        logger.info("Parsing chromossome size file '%s' and creating dic." % chrmFile)
    try:
        with open(chrmFile, 'r') as f:
            chrmSizes = {}
            for line in enumerate(f):
                row = line[1].strip().split('\t')
                chrmSizes[str(row[0])] = int(row[1])
        return chrmSizes
    except IOError:
        if not s:
            logger.error(" '%s' file not openable or doesn't exist." % chrmFile)
            print("%s file not openable or doesn't exist." % chrmFile)
        sys.exit(0)

def parseGFF(infile, chrmSizes, promSize, operons, operonDist):
    """
    Parses GFF file and annotates genomic features.
    Promoters 
    Outputs list of lists with annotation for each feature.
    """
    if not s:
        logger.info("Parsing gff file '%s' and creating annotation." % infile)
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

                    # check if same chromossome as before
                    if (prev_chrm == chrm) and (prev_chrm != ""):
                        if trigger:
                            if strand == "+":
                                # two genes head to head situation
                                # test if there's space for both promoters
                                if prev_end + 1 + (promSize * 2) + 1 < start:
                                    # TSS of last gene
                                    cur_loc = "TSS"
                                    output.append([prev_chrm, prev_end, prev_end + 1, cur_loc, prev_gene])
                                    # promoter of last gene
                                    istart = prev_end + 1
                                    iend = prev_end + 1 + promSize
                                    cur_loc = "Promoter"
                                    output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                    # intergenic space
                                    cur_loc = "Intergenic"
                                    istart = prev_end + 1 + promSize + 1
                                    iend = start - 1 - promSize - 1
                                    output.append([chrm, istart, iend, cur_loc, "."])
                                    # New gene promoter
                                    istart = start - 1 - promSize - 1
                                    iend = start - 1
                                    cur_loc = "Promoter"
                                    output.append([chrm, istart, iend, cur_loc, gene])
                                    # TSS of current gene
                                    cur_loc = "TSS"
                                    output.append([prev_chrm, start, start + 1, cur_loc, gene])
                                    
                                else:
                                    # there's not enough space for both promoters, allow overlapping promoters
                                    # Another option: divide space equally, non-overlapping
                                    # TSS of last gene
                                    cur_loc = "TSS"
                                    output.append([prev_chrm, prev_end, prev_end + 1, cur_loc, prev_gene])
                                    # promoter of last gene
                                    istart = prev_end + 1
                                    iend = prev_end + 1 + promSize
                                    cur_loc = "Promoter"
                                    output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                    # intergenic space
                                    cur_loc = "Intergenic"
                                    istart = prev_end + 1 + promSize + 1
                                    iend = start - 1 - promSize - 1
                                    output.append([chrm, istart, iend, cur_loc, "."])
                                    # New gene promoter
                                    istart = start - 1 - promSize - 1
                                    iend = start - 1
                                    cur_loc = "Promoter"
                                    output.append([chrm, istart, iend, cur_loc, gene])
                                    # TSS of current gene
                                    cur_loc = "TSS"
                                    output.append([prev_chrm, start, start + 1, cur_loc, gene])
                                    
                                trigger = False
                            else:
                                # two genes in neg position

                                if operons:
                                    # test if in operon
                                    if prev_end + operonDist > start:
                                        # not operon
                                        # test if there's space for promoter of the last gene
                                        if prev_end + 1 + promSize < start:
                                            # there's enough space
                                            # TSS of last gene
                                            cur_loc = "TSS"
                                            output.append([prev_chrm, prev_end - 1, prev_end, cur_loc, prev_gene])
                                            # promoter of last gene
                                            istart = prev_end + 1
                                            iend = prev_end + 1 + promSize
                                            cur_loc = "Promoter"
                                            output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                            # add intergenic space following
                                            cur_loc = "Intergenic"
                                            istart = prev_end + 1 + promSize + 1
                                            iend = start - 1
                                            output.append([prev_chrm, istart, iend, cur_loc, "."])
                                                                                  
                                        else:
                                            # there's not enough space, add remaining until next gene
                                            # TSS of last gene
                                            cur_loc = "TSS"
                                            output.append([prev_chrm, prev_end - 1, prev_end, cur_loc, prev_gene])
                                            # promoter of last gene
                                            istart = prev_end + 1
                                            iend = start - 1
                                            cur_loc = "Promoter"
                                            output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                            # skip intergenic space
                                    else:
                                        # in operon
                                        # last promoter skipped
                                        # last TSS skipped
                                        # fill with intergenic space
                                        cur_loc = "Intergenic"
                                        istart = prev_end + 1
                                        iend = start - 1
                                        output.append([prev_chrm, istart, iend, cur_loc, "."])
                                else:
                                    # test if there's space for promoter of the last gene
                                    if prev_end + 1 + promSize < start:
                                        # there's enough space
                                        # TSS of last gene
                                        cur_loc = "TSS"
                                        output.append([prev_chrm, prev_end - 1, prev_end, cur_loc, prev_gene])
                                        # promoter of last gene
                                        istart = prev_end + 1
                                        iend = prev_end + 1 + promSize
                                        cur_loc = "Promoter"
                                        output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                        # add intergenic space following
                                        cur_loc = "Intergenic"
                                        istart = prev_end + 1 + promSize + 1
                                        iend = start - 1
                                        output.append([prev_chrm, istart, iend, cur_loc, "."])
                                                                              
                                    else:
                                        # there's not enough space, add remaining until next gene
                                        # TSS of last gene
                                        cur_loc = "TSS"
                                        output.append([prev_chrm, prev_end - 1, prev_end, cur_loc, prev_gene])
                                        # promoter of last gene
                                        istart = prev_end + 1
                                        iend = start - 1
                                        cur_loc = "Promoter"
                                        output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                        # skip intergenic space

                                trigger = True
                        else:
                            # no trigger
                            if strand == "+":
                                # two positive genes situation
                                
                                if operons:
                                    # test if in operon
                                    if prev_end + operonDist > start:
                                        # not operon
                                        # test if there's space for promoter
                                        if start - 1 - promSize > prev_end:
                                            # add intergenic
                                            cur_loc = "Intergenic"
                                            istart = prev_end + 1
                                            iend = start - 1 - promSize
                                            output.append([chrm, istart, iend, cur_loc, "."])
                                            # promoter of last
                                            istart = start - 1 - promSize
                                            iend = start - 1
                                            cur_loc = "Promoter"
                                            output.append([chrm, istart, iend, cur_loc, gene])
                                            # TSS of current gene
                                            cur_loc = "TSS"
                                            output.append([prev_chrm, start, start + 1, cur_loc, gene])
                                            
                                        else:
                                            # if there's not enough space, fill space with promoter
                                            # skip intergenic space
                                            istart = prev_end + 1
                                            iend = start - 1
                                            cur_loc = "Promoter"
                                            output.append([chrm, istart, iend, cur_loc, gene])
                                            # TSS of current gene
                                            cur_loc = "TSS"
                                            output.append([prev_chrm, start, start + 1, cur_loc, gene])

                                else:
                                    # in operon
                                        # fill with intergenic space 
                                        # skip promoter
                                        # skip TSS
                                        cur_loc = "Intergenic"
                                        istart = prev_end + 1
                                        iend = start - 1
                                        output.append([prev_chrm, istart, iend, cur_loc, "."])
                                    
                                trigger = False
                            else:
                                # two gene endings situation
                                # fill with intergenic space
                                cur_loc = "Intergenic"
                                istart = prev_end + 1
                                iend = start - 1
                                output.append([prev_chrm, istart, iend, cur_loc, "."])
                                
                                trigger = True                                  
                    else:
                        # chromossome change                        
                        if prev_chrm != "": # avoid inititallized chr
                            if trigger:
                                # check if there's enough space for promoter of the last chr
                                if prev_end + 1 + promSize < chrmSizes[prev_chrm]:
                                    # TSS of last gene
                                    cur_loc = "TSS"
                                    output.append([prev_chrm, prev_end, prev_end + 1, cur_loc, prev_gene])
                                    # promoter of last gene
                                    istart = prev_end + 1
                                    iend = prev_end + 1 + promSize
                                    cur_loc = "Promoter"
                                    output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                                    # add intergenic until end of chr
                                    cur_loc = "Intergenic"
                                    istart = prev_end + 1 + promSize
                                    iend = chrmSizes[prev_chrm]
                                    output.append([prev_chrm, istart, iend, cur_loc, "."])
                                else:
                                    # if not fill rest of chromossome with promoter
                                    # TSS of last gene
                                    cur_loc = "TSS"
                                    output.append([prev_chrm, prev_end - 1, prev_end, cur_loc, prev_gene])
                                    # promoter of last gene
                                    if prev_end != chrmSizes[prev_chrm]:
                                        # avoid going beyond genome boundaries
                                        istart = prev_end + 1
                                        iend = chrmSizes[prev_chrm]
                                        cur_loc = "Promoter"
                                        output.append([prev_chrm, istart, iend, cur_loc, prev_gene])
                            else:
                                # fill rest of last chromossome with intergenic
                                cur_loc = "Intergenic"
                                istart = prev_end + 1
                                iend = chrmSizes[prev_chrm]
                                output.append([prev_chrm, istart, iend, cur_loc, "."])
                        
                        # operate on new chromossome
                        if strand == "+":
                            # test if there's space for promoter
                            if start - 1 - promSize > 1:
                                # there's enough space
                                # add intergenic space from beggining of chr
                                cur_loc = "Intergenic"
                                istart = 1
                                iend = start - 1 - promSize
                                output.append([chrm, istart, iend, cur_loc, "."])
                                # add promoter
                                if start > 2:
                                    # avoid going beyond genome boundaries
                                    istart = start - 1 - promSize
                                    iend = start - 1
                                    cur_loc = "Promoter"
                                    output.append([chrm, istart, iend, cur_loc, gene])
                                # TSS of current gene
                                cur_loc = "TSS"
                                output.append([chrm, start, start + 1, cur_loc, gene])
                                
                            else:
                                # skip intergenic
                                if start > 2:
                                    # avoid going beyond genome boundaries
                                    istart = 1
                                    iend = start - 1 # make sure it's > 1
                                    cur_loc = "Promoter"
                                    output.append([chrm, istart, iend, cur_loc, gene])
                                # TSS of current gene
                                cur_loc = "TSS"
                                output.append([chrm, start, start + 1, cur_loc, gene])
                                
                            trigger = False
                        else:
                            # add intergenic space from beggining of chr
                            cur_loc = "Intergenic"
                            istart = 1
                            iend = start - 1
                            output.append([chrm, istart, iend, cur_loc, "."])
                            trigger = True
                            
                        
                    start = istart
                    end = iend

                elif cur_line == "mRNA":
                    pass

                elif cur_line == "CDS":
                    # annotate intron if two consecutive CDSs or UTR preceding CDS not consecutive
                    if (prev_loc == "CDS" and prev_end + 1 != start) or ("UTR" in prev_loc and prev_end + 1 != start):
                        
                        cur_loc = "Intron"
                        istart = prev_end + 1
                        # make sure it has length >= 1
                        if start - 1 > istart:
                            iend = start - 1
                        else:
                            iend = istart
                        output.append([chrm, istart, iend, cur_loc, gene])

                    # annotate CDS
                    cur_loc = "CDS"
                    output.append([chrm, start, end, cur_loc, gene])

                elif cur_line == "UTR":
                    # annotate intron if two consecutive UTRs or CDS precedding UTR non-overlapping
                    if (prev_loc == "CDS" and prev_end + 1 != start) or ("UTR" in prev_loc and prev_end + 1 != start):
                        cur_loc = "Intron"
                        istart = prev_end + 1
                        iend = start - 1
                        output.append([chrm, istart, iend, cur_loc, gene])

                    # annotate 5' UTR
                    if (strand == "+" and prev_loc == "TSS") or (strand == "-" and (prev_loc == "CDS" or "UTR" in prev_loc)):
                        cur_loc = "5'UTR"
                        output.append([chrm, start, end, cur_loc, gene])

                    # annotate 3'
                    elif (strand == "-" and (prev_loc == "Intergenic" or prev_loc == "3'UTR")) or (strand == "+" and (prev_loc == "CDS")):
                        cur_loc = "3'UTR"
                        output.append([chrm, start, end, cur_loc, gene])

                prev_line = cur_line
                prev_chrm = chrm
                prev_end = end
                prev_loc = cur_loc
                prev_gene = gene

        # Annotate chrsm without genes
        allChrms = set(chrmSizes.keys())
        curChrs = set([c[0] for c in output])

        for c in allChrms.difference():
            output.append([c, 1, chrmSizes[c], "Intergenic"])

        return output#.sort()

    except IOError:
        if not s:
            logger.error(" '%s' file not openable or doesn't exist." % infile)
            print(" '%s' file not openable or doesn't exist." % infile)
        sys.exit(0)

def writeOutput(output, outFile):
    """
    Writes output (list) tab-delimited
    to outFile or to stout if outFile is not specified.
    """
    try:
        if outFile != "":
            if not s:
                logger.info("Writing output to '%s'." % outFile)
            # write to file
            with open(outFile, 'wb') as f:
                wr = csv.writer(f, delimiter='\t')
                for line in output:
                    wr.writerow(line)
        else:
            if not s:
                logger.info("Writing output to stdout.")
            wr = csv.writer(sys.stdout, delimiter='\t')
            for line in output:
                wr.writerow(line)
    except IOError:
        if not s:
            logger.error(" '%s' file not writable." % outFile)
            print(" '%s' file not writable." % outFile)
        sys.exit(0)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)