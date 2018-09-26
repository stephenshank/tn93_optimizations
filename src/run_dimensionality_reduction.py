import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA


def dimensionality_reduction(input, output, dimension, method):
  dimension = int(dimension)
  kmer_df = pd.read_csv(input, index_col=0)
  scaled_kmer = StandardScaler().fit_transform(kmer_df)
  if method == 'tsne':
    embedded_kmer = TSNE(n_components=dimension, n_iter=10000, random_state=0).fit_transform(scaled_kmer)
  else:
    embedded_kmer = PCA(n_components=dimension, random_state=0).fit_transform(scaled_kmer)
  results = {'header': kmer_df.index}
  for i in range(dimension):
      results['x'+str(i)] = embedded_kmer[:, i]
  pd.DataFrame(results).to_csv(output, index=False)

