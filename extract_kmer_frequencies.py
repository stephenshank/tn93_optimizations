import argparse
import itertools as it

import numpy as np
from Bio import SeqIO


parser = argparse.ArgumentParser(
  description="Extract k-mer counts (or frequencies) from a fasta file."
)

parser.add_argument(
  '-i', '--input',
  metavar='INPUT',
  type=str,
  help='Input fasta file.',
  required=True,
)

parser.add_argument(
  '-o', '--output',
  metavar='OUTPUT',
  type=str,
  help='Output CSV file.',
  required=True,
)

parser.add_argument(
  '-k', '--kmer',
  metavar='KMER',
  type=str,
  help='Length of k-mer.',
  required=True,
)

parser.add_argument(
  '-H', '--header',
  help='Whether or not header should be written.',
  default=True,
  action='store_true'
)

args = parser.parse_args()

k = int(args.kmer)
N = 4**k
kmer_tuples = it.product('ACGT', repeat=k)
kmer_to_index = { ''.join(kmer): index for index, kmer in enumerate(kmer_tuples) }
output_file = open(args.output, 'w')

if args.header:
  kmers = [''.join(kmer) for kmer in it.product('ACGT', repeat=k)]
  output_file.write(','.join(['header']+kmers))

records = SeqIO.parse(args.input, 'fasta')
bins = np.arange(N+1)
for record in records:
  counts = np.zeros(N, dtype=np.int)
  ungapped_string = str(record.seq.ungap('-'))
  current_kmers = [ungapped_string[i:i+k] for i in range(len(ungapped_string)-k+1)]
  indices = [kmer_to_index[kmer] for kmer in current_kmers if kmer in kmer_to_index]
  kmer_counts, _ = np.histogram(indices, bins=bins)
  output_file.write('\n'+record.name+',')
  output_file.write(','.join([str(i) for i in kmer_counts]))

