import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Parsing arguments
parser = argparse.ArgumentParser(description='Remove gap in "fasta" file.')
parser.add_argument('--input', '-i', required=True, help='Input file')
parser.add_argument('--output', '-o', required=True, help='Output file')
d = vars(parser.parse_args())
input_file, output_file = d['input'], d['output']

with open(output_file, "w") as output:
    for record in SeqIO.parse(input_file, "fasta"):
        record.seq = Seq(str(record.seq).replace("-", ""))
        SeqIO.write(record, output, "fasta")
