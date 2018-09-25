import argparse

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
  help='Output fasta file.',
  required=True,
)

parser.add_argument(
  '-n', '--number',
  metavar='NUMBER',
  type=str,
  help='Number of entries to subset.',
  required=True,
)

args = parser.parse_args()

records = list(SeqIO.parse(args.input, 'fasta'))
np.random.seed(0)
indices = np.random.choice(len(records), int(args.number), replace=False)
SeqIO.write(
  [record for i, record in enumerate(records) if i in indices],
  args.output,
  'fasta'
)

