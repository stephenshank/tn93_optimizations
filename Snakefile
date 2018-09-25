from src import *

rule subset:
  input:
    fasta="input/HIV-LANL-unaligned.fasta"
  output:
    fasta="output/HIV-LANL-unaligned_subset-{subset}.fasta"
  run:
    subset(input.fasta, output.fasta, wildcards.subset)

rule alignment:
  input:
    rules.subset.output.fasta
  output:
    "output/HIV-LANL_subset-{subset}.fasta"
  shell:
    "mafft {input} > {output}"

rule tn93_distances:
  input:
    rules.alignment.output
  output:
    "output/HIV-LANL_subset-{subset}_tn93.csv"
  shell:
    "tn93 -o {output} -t 1 {input}"

rule kmers:
  input:
    fasta=rules.subset.output.fasta
  output:
    csv="output/HIV-LANL_subset-{subset}_{k}-mers.csv"
  run:
    extract_kmer_frequencies(input.fasta, output.csv, wildcards.k)

rule tsne_distances:
  input:
    csv=rules.kmers.output.csv
  output:
    csv="output/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim_tsne.csv"
  run:
    tsne_distances(input.csv, output.csv, wildcards.dimension)

rule pairwise_information:
  input:
    rules.tn93_distances.output,
    rules.tsne_distances.output.csv
  output:
    "output/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.csv"
  run:
    pairwise_information(wildcards.subset, wildcards.k, wildcards.dimension)

rule plot:
  input:
    rules.pairwise_information.output
  output:
    "plots/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.png"
  run:
    plot(wildcards.subset, wildcards.k, wildcards.dimension)

SUBSETS = ["1000"]
KS = ["3", "4", "5"]
DIMENSIONS = ["2", "3"]
rule all:
  input:
    expand("plots/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.png", subset=SUBSETS, k=KS, dimension=DIMENSIONS)
