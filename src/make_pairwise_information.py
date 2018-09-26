import csv

import numpy as np
import pandas as pd


def pairwise_information(tn93_input, dr_input, output, dimension):
  dimension = int(dimension)
  reduced_dict = {}
  with open(dr_input) as csv_file:
    reduced_csv = csv.DictReader(csv_file)
    for row in reduced_csv:
      reduced_dict[row['header']] = {'x'+str(i): float(row['x'+str(i)]) for i in range(dimension)}

  tn93_df = pd.read_csv(tn93_input).rename(columns={'Distance': 'TN93 Distance'})
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

  tn93_df.to_csv(output, index=False)

