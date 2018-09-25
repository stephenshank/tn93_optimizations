import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def plot(subset, k, char, dimension):
  parameters = ( subset, k, char, dimension )

  pairwise_filename = "output/HIV-LANL_subset-%s_%s-mers_%s-char_%s-dim.csv" % parameters
  tn93_df = pd.read_csv(pairwise_filename)
  plot_filename = "plots/HIV-LANL_subset-%s_%s-mers_%s-char_%s-dim.png" % parameters
  sns.jointplot(x="TN93 Distance", y="TSNE Euclidean Distance", data=tn93_df, kind='hex')
  plt.savefig(plot_filename)

