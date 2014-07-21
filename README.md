annotateFromGFF.py
=======

Creates functional annotation of complete genome based on features on GFF file.

Reports coding-sequence (CDS) and untranslated-regions (UTR) and extracts features such as introns, intergenic space and promoters. It also distinguishes between 5' and 3' UTRs.

Promoters and intergenic space are annotated dynamically based on a specified promoter size but rezisable to fit genome boundaries and small intergenic space.

# Usage

`python annotateFromGFF.py [OPTIONS] file.gff chrmSizes.tsv > annotation.bed`

## Positional arguments (required)

gff - GFF file with annotation.

chrmFile - Tab-delimited file with sizes of each chromossome (chr:size).

## Optional arguments

`-o`, `--outfile` - Specifies the name of the output file. If not specified, will output to stdout.

`-p`, `--promoterSize` - Average size of promoter elements. Dynamically resizable.

`-v` - Print debugging information

# Promoter annotation
Promoters are defined as a region upstream of transcription start sites (TSSs) by a fixed length (-p argument) but are resized to less if that distance would overlap a chromossome boundary or other gene.

Intergenic regions are resized accordingly to the complementary space between two genes minus the promoter space.

