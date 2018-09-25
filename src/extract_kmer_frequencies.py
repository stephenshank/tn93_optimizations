import itertools as it

import numpy as np
from Bio import SeqIO


def extract_kmer_frequencies(input, output, k):
  k = int(k)
  N = 4**k
  kmer_tuples = it.product('ACGT', repeat=k)
  kmer_to_index = { ''.join(kmer): index for index, kmer in enumerate(kmer_tuples) }
  output_file = open(output, 'w')
  kmers = [''.join(kmer) for kmer in it.product('ACGT', repeat=k)]
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

