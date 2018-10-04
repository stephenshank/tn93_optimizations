import itertools as it
import csv

import pandas as pd
import numpy as np
from Bio import SeqIO


def extract_kmer_frequencies(input, output, k, characters='nuc'):
  k = int(k)
  if characters == 'nuc':
    characters = 'AGCT'
  else:
    characters = 'ACDEFGHIKLMNPQRSTVWY'
  N = len(characters)**k
  kmer_tuples = it.product(characters, repeat=k)
  kmer_to_index = { ''.join(kmer): index for index, kmer in enumerate(kmer_tuples) }
  output_file = open(output, 'w')
  kmers = [''.join(kmer) for kmer in it.product(characters, repeat=k)]
  output_file.write(','.join(['header']+kmers))

  records = SeqIO.parse(input, 'fasta')
  bins = np.arange(N+1)
  for record in records:
    counts = np.zeros(N, dtype=np.int)
    ungapped_string = str(record.seq.ungap('-'))
    current_kmers = [ungapped_string[i:i+k] for i in range(len(ungapped_string)-k+1)]
    indices = [kmer_to_index[kmer] for kmer in current_kmers if kmer in kmer_to_index]
    kmer_counts, _ = np.histogram(indices, bins=bins)
    output_file.write('\n'+record.name+',')
    output_file.write(','.join([str(i) for i in kmer_counts]))


def window(input, output, length):
  with open(output, 'w') as csv_file:
    length = int(length)
    if length == 0:
      return
    records = SeqIO.parse(input, 'fasta')
    reference = next(records)
    windows = np.arange(0, len(reference), length)
    refnp = np.array(list(str(reference.seq)), dtype='<U1')
    writer = csv.writer(csv_file)
    writer.writerow(['header']+['window%d' % i for i in range(len(windows)-1)])
    for record in records:
      recordnp = np.array(list(str(record.seq)), dtype='<U1')
      difference_cumsum = np.cumsum(refnp != recordnp)
      window_differences = difference_cumsum[windows[1:]]-difference_cumsum[windows[0:-1]]
      writer.writerow([record.name]+list(window_differences))


def concat(kmer_input, window_input, output):
  length = int(output.split('_')[-1].split('.')[0])
  kmers = pd.read_csv(kmer_input, index_col=0)
  if length == 0:
    result = kmers
  else:
    window = pd.read_csv(window_input, index_col=0)
    result = pd.concat([kmers, window], axis=1)
  result.to_csv(output)

