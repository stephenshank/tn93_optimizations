import argparse

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(
  description="Extract k-mer counts (or frequencies) from a fasta file."
)

parser.add_argument(
  '-k', '--kmer',
  metavar='KMER',
  type=int,
  help='K to determine kmers.',
  required=True,
)

parser.add_argument(
  '-s', '--subset',
  metavar='SUBSET',
  type=int,
  help='Amount of subsetted sequences.',
  required=True,
)

parser.add_argument(
  '-d', '--dimension',
  metavar='DIMENSION',
  type=int,
  help='Number of reduced dimension.',
  required=True,
)

args = parser.parse_args()
parameters = ( args.subset, args.kmer, args.dimension )

pairwise_filename = "output/HIV-LANL_subset-%d_%d-mers_%d-dim.csv" % parameters
tn93_df = pd.read_csv(pairwise_filename)
plot_filename = "plots/HIV-LANL_subset-%d_%d-mers_%d-dim.png" % parameters
sns.jointplot(x="TN93 Distance", y="TSNE Euclidean Distance", data=tn93_df, kind='hex')
plt.savefig(plot_filename)

