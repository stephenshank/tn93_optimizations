import argparse

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE

parser = argparse.ArgumentParser(
  description="Run t-SNE dimensionality reduction on a CSV."
)

parser.add_argument(
  '-i', '--input',
  metavar='INPUT',
  type=str,
  help='Input CSV file.',
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
  '-d', '--dimension',
  metavar='DIMENSION',
  type=str,
  help='Number of reduced dimensions.',
  required=True,
)

args = parser.parse_args()
dimension = int(args.dimension)
kmer_df = pd.read_csv(args.input, index_col=0)
scaled_kmer = StandardScaler().fit_transform(kmer_df)
embedded_kmer = TSNE(n_components=dimension, n_iter=10000, random_state=0).fit_transform(scaled_kmer)
results = {'header': kmer_df.index}
for i in range(dimension):
    results['x'+str(i)] = embedded_kmer[:, i]
pd.DataFrame(results).to_csv(args.output, index=False)

