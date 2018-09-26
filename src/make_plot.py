import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def plot(input, output):
  tn93_df = pd.read_csv(input)
  max_tsne_distance = tn93_df['TSNE Euclidean Distance'].max()
  tn93_df['Normalized TSNE Euclidean Distance'] = tn93_df['TSNE Euclidean Distance'] / max_tsne_distance
  below_threshold = tn93_df['TN93 Distance'] <= 1
  ax = sns.jointplot(
    x="TN93 Distance",
    y="Normalized TSNE Euclidean Distance",
    data=tn93_df[below_threshold],
    kind='hex'
  )
  plt.subplots_adjust(top=0.9)
  _, subset, gene, k, info = input.split('/')
  _, char, dimension, method = info.split('_')
  method = method.split('.')[0]
  title_params = ( gene, subset, char, k, method )
  title_str = '%s - %s %s %s-mers, method=%s, reduced' % title_params
  ax.fig.suptitle(title_str)
  plt.savefig(output)

