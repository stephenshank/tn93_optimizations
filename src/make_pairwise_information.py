import csv

import numpy as np
import pandas as pd


def pairwise_information(subset, k, char, dimension):
  subset = int(subset)
  k = int(k)
  dimension = int(dimension)
  parameters = ( subset, k, char, dimension )

  reduced_csv_filename = "output/HIV-LANL_subset-%d_%d-mers_%s-char_%d-dim_tsne.csv" % parameters
  reduced_dict = {}
  with open(reduced_csv_filename) as csv_file:
    reduced_csv = csv.DictReader(csv_file)
    for row in reduced_csv:
      reduced_dict[row['header']] = {'x'+str(i): float(row['x'+str(i)]) for i in range(dimension)}

  tn93_distance_filename = "output/HIV-LANL_subset-%d_tn93.csv" % parameters[0]
  tn93_df = pd.read_csv(tn93_distance_filename).rename(columns={'Distance': 'TN93 Distance'})
  for i in range(dimension):
      tn93_df['x'+str(i)] = np.zeros(len(tn93_df))
      tn93_df['y'+str(i)] = np.zeros(len(tn93_df))
  for key, val in reduced_dict.items():
    for i in range(dimension):
        tn93_df.loc[tn93_df.ID1 == key, 'x'+str(i)] = val['x'+str(i)]
        tn93_df.loc[tn93_df.ID2 == key, 'y'+str(i)] = val['x'+str(i)]

  tn93_df['TSNE Euclidean Distance'] = np.zeros(len(tn93_df))
  for i in range(dimension):
      tn93_df['TSNE Euclidean Distance'] += (tn93_df['x'+str(i)]-tn93_df['y'+str(i)])**2
  tn93_df['TSNE Euclidean Distance'] = np.sqrt(tn93_df['TSNE Euclidean Distance'])

  pairwise_filename = "output/HIV-LANL_subset-%d_%d-mers_%s-char_%d-dim.csv" % parameters
  tn93_df.to_csv(pairwise_filename, index=False)
