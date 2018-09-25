from src import *

rule subset:
  input:
    fasta="input/HIV-LANL-unaligned.fasta"
  output:
    fasta="output/HIV-LANL-unaligned_subset-{subset}.fasta"
  run:
    subset(input.fasta, output.fasta, wildcards.subset)

rule alignment_bam:
  input:
    rules.subset.output.fasta
  output:
    "output/HIV-LANL_subset-{subset}.bam"
  shell:
    "bealign -r HXB2_env {input} {output}"

rule alignment:
  input:
    rules.alignment_bam.output
  output:
    "output/HIV-LANL_subset-{subset}_nuc-char.fasta"
  shell:
    "bam2msa {input} {output}"

rule protein_alignment:
  input:
    rules.alignment.output
  output:
    "output/HIV-LANL_subset-{subset}_aa-char.fasta"
  shell:
    "translate {input} > {output}"

rule tn93_distances:
  input:
    rules.alignment.output
  output:
    "output/HIV-LANL_subset-{subset}_tn93.csv"
  shell:
    "tn93 -o {output} -t 1 {input}"

rule kmers:
  input:
    fasta="output/HIV-LANL_subset-{subset}_{char}-char.fasta"
  output:
    csv="output/HIV-LANL_subset-{subset}_{k}-mers_{char}-char.csv"
  run:
    extract_kmer_frequencies(input.fasta, output.csv, wildcards.k)

rule tsne_distances:
  input:
    csv=rules.kmers.output.csv
  output:
    csv="output/HIV-LANL_subset-{subset}_{k}-mers_{char}-char_{dimension}-dim_tsne.csv"
  run:
    tsne_distances(input.csv, output.csv, wildcards.dimension)

rule pairwise_information:
  input:
    rules.tn93_distances.output,
    rules.tsne_distances.output.csv
  output:
    "output/HIV-LANL_subset-{subset}_{k}-mers_{char}-char_{dimension}-dim.csv"
  run:
    pairwise_information(wildcards.subset, wildcards.k,  wildcards.char, wildcards.dimension)

rule plot:
  input:
    rules.pairwise_information.output
  output:
    "plots/HIV-LANL_subset-{subset}_{k}-mers_{char}-char_{dimension}-dim.png"
  run:
    plot(wildcards.subset, wildcards.k, wildcards.char, wildcards.dimension)

SUBSETS = ["100"]
NUCLEOTIDE_KS = ["3", "4", "5"]
DIMENSIONS = ["2", "3"]
rule all:
  input:
    expand("plots/HIV-LANL_subset-{subset}_{k}-mers_nuc-char_{dimension}-dim.png", subset=SUBSETS, k=NUCLEOTIDE_KS, dimension=DIMENSIONS),
    expand("plots/HIV-LANL_subset-{subset}_2-mers_aa-char_{dimension}-dim.png", subset=SUBSETS, dimension=DIMENSIONS),
