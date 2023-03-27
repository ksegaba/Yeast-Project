"""
Take each nucleotide sequences from a fasta file and translate it to protein
and generate a new fasta file.

Inputs:
- CDS fasta file

Outputs:
- Protein fasta file
"""

import argparse, sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def make_protein_record(nucl_record):
    # Translates nucleotide sequence to protein sequence
    return SeqRecord(
        seq=nucl_record.seq.translate(cds=True),
        id=nucl_record.id,
        description='translation from CDS, using default table')

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
        description='Translate CDS fasta file to protein')
    
    # Required input
    req_group = parser.add_argument_group(title='REQUIRED INPUT')
    req_group.add_argument('-cds', help='CDS fasta file', required=True)
    req_group.add_argument(
        '-save', help='Save name of protein fasta file', required=True)

    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()  # Read arguments

    # Translate CDS fasta to protein
    proteins = (
        make_protein_record(nucl_record)
        for nucl_record in SeqIO.parse(args.cds, 'fasta'))

    # Write new protein fasta file
    SeqIO.write(proteins, args.save, 'fasta')

if __name__ == '__main__':
    main()
