#!/usr/bin/env python

# Slice subregions out of a genbank between 2 genes.
# Usage:
#   slice_genes.py -i infile.gbk -r gene1:gene2 -s strain -o outfile
# From: https://www.biostars.org/p/340270/

from Bio import SeqIO
import sys, argparse, os

def get_args():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description='Subset genbanks between 2 genes.')
        parser.add_argument('-i', '--infile', action='store',
                            help='Input genbank file to slice.')
        parser.add_argument('-r', '--range', action='store',
                            help='The 2 gene LOCUS TAGS to slice between')
        parser.add_argument('-o', '--outfile', action='store',
                            help='Output file name prefix.')
        parser.add_argument('-s', '--strain', action='store',
                            help='Strain name for gene prefix.')
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)
    except:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")

    return parser.parse_args()


def main():
    args = get_args()
    genes = args.range.split(':')
    sys.stderr.write('Fetching subsequence betweeen {} and {}.\n'.format(genes[0], genes[1]))

    records = SeqIO.parse(args.infile, 'genbank')
    buffer = "Name,Feature\n"
    for record in records:
        sys.stderr.write('Operating on {}.\n'.format(record.id))
        loci = [feat for feat in record.features if feat.type == "CDS"]
        try:
            start = min([int(l.location.start) for l in loci if l.qualifiers['locus_tag'][0] in genes])
            end = max([int(l.location.end) for l in loci if l.qualifiers['locus_tag'][0] in genes])
        except ValueError:
            sys.stderr.write('No indices returned for those loci, assume they don\'t feature in this record, moving on...\n.')
            continue
        try:
            (start and end)
            subrecord = record[start:end]
            # subrecord.molecule_type = "DNA"
            # print(vars(subrecord))
            seq = subrecord.seq
            subrecord.id = os.path.basename(args.outfile)
            subrecord.description = ""
            SeqIO.write(subrecord, args.outfile + ".fasta", "fasta")
            for feature in subrecord.features:
                if feature.type == "gene":
                    loc_str = str(feature.location)
                    f_start, f_stop = loc_str[1:].split("]")[0].split(":")
                    dir = loc_str[-2]
                    f_seq = seq[int(f_start):int(f_stop)]
                    if dir == "-":
                        f_seq = f_seq.reverse_complement()
                    buffer = F"{buffer}{args.strain}_{feature.qualifiers['locus_tag'][0]},{f_seq}\n"
                    print(F"Start: {f_start}, Stop: {f_stop}, Strand: {dir}, Type: {feature.type}, Name: {feature.qualifiers['locus_tag'][0]}")
            with open(args.outfile + ".csv", "w") as ofile:
                ofile.write(buffer)
            break
        except NameError:
            sys.stderr.write('Didn\'t get any indices even though the genes seemed to match. Couldn\'t slice.\n')
            pass


if __name__ == '__main__':
    main()
