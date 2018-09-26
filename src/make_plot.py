import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def plot(subset, k, char, dimension):
  parameters = ( subset, k, char, dimension )

  pairwise_filename = "output/%s/%s/pairwise_%s_%sd.csv" % parameters
  tn93_df = pd.read_csv(pairwise_filename)
  max_tsne_distance = tn93_df['TSNE Euclidean Distance'].max()
  tn93_df['Normalized TSNE Euclidean Distance'] = tn93_df['TSNE Euclidean Distance'] / max_tsne_distance
  plot_filename = "output/%s/%s/plots/%s_%sd.png" % parameters
  below_threshold = tn93_df['TN93 Distance'] <= 1
  ax = sns.jointplot(
    x="TN93 Distance",
    y="Normalized TSNE Euclidean Distance",
    data=tn93_df[below_threshold],
    kind='hex'
  )
  plt.subplots_adjust(top=0.9)
  title_params = ( subset, char, k )
  title_str = '%s %s %s-mers, method=t-sne, reduced' % title_params
  ax.fig.suptitle(title_str)
  plt.savefig(plot_filename)

