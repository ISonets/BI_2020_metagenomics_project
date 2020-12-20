import sys
from Bio import SeqIO
from Bio.Seq import Seq

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, "w") as output:
    for record in SeqIO.parse(input_file, "fasta"):
        record.seq = Seq(str(record.seq).replace("-", ""))
        SeqIO.write(record, output, "fasta")