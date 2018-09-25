import numpy as np
from Bio import SeqIO


def subset(input, output, number):
    records = list(SeqIO.parse(input, 'fasta'))
    np.random.seed(0)
    indices = np.random.choice(len(records), int(number), replace=False)
    SeqIO.write(
      [record for i, record in enumerate(records) if i in indices],
      output,
      'fasta'
    )

