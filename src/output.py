import csv

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def plot(input, output, distance):
  tn93_df = pd.read_csv(input)
  max_tsne_distance = tn93_df['TSNE Euclidean Distance'].max()
  tn93_df['Normalized TSNE Euclidean Distance'] = tn93_df['TSNE Euclidean Distance'] / max_tsne_distance
  max_tsne_distance = tn93_df['TSNE L1 Distance'].max()
  tn93_df['Normalized TSNE L1 Distance'] = tn93_df['TSNE L1 Distance'] / max_tsne_distance
  below_threshold = tn93_df['TN93 Distance'] <= .015
  ax = sns.jointplot(
    x="TN93 Distance",
    y="Normalized TSNE %s Distance" % distance,
    data=tn93_df[below_threshold],
    kind='hex'
  )
  plt.subplots_adjust(top=0.9)
  _, subset, gene, k, info = input.split('/')
  _, char, dimension, method = info.split('_')
  method = method.split('.')[0]
  title_params = ( gene, subset, char, k, method, distance )
  title_str = '%s - %s %s %s-mers, method=%s, %s' % title_params
  ax.fig.suptitle(title_str)
  plt.savefig(output)

def table(output, subsets):
  _, gene, k, info = output.split('/')
  char, dimension, method = info.split('_')
  method = method.split('.')[0]
  with open(output, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['subset', 'distance', 'max', '99th p', '95th p', '75th p'])
    for subset in subsets:
      filename_params = ( subset, gene, k, char, dimension, method )
      filename = 'output/%d/%s/%s/pairwise_%s_%s_%s.csv' % filename_params
      df = pd.read_csv(filename)
      total_pairs = len(df)
      for distance in ['Euclidean', 'L1']:
        distance_string = 'TSNE %s Distance' % distance
        rd_desired = df[df['TN93 Distance'] < .015][distance_string]
        dp_max = rd_desired.max()
        max_pairs = (df[distance_string] < dp_max).sum()/total_pairs
        dp_99 = rd_desired.quantile(.99)
        dp_99_pairs = (df[distance_string] < dp_99).sum()/total_pairs
        dp_95 = rd_desired.quantile(.95)
        dp_95_pairs = (df[distance_string] < dp_95).sum()/total_pairs
        dp_75 = rd_desired.quantile(.75)
        dp_75_pairs = (df[distance_string] < dp_75).sum()/total_pairs
        row = [
          subset,
          distance,
          max_pairs,
          dp_99_pairs,
          dp_95_pairs,
          dp_75_pairs
        ]
        writer.writerow(row)

