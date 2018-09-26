from src import *

rule uncompressed_data:
  input:
    "input/HIV-LANL-unaligned.fasta.tar.gz"
  output:
    "input/HIV-LANL-unaligned.fasta"
  shell:
    "tar xvzf input/HIV-LANL-unaligned.fasta.tar.gz"

rule subset:
  input:
    fasta=rules.uncompressed_data.output[0]
  output:
    fasta="output/{subset}/HIV.fasta"
  run:
    subset(input.fasta, output.fasta, wildcards.subset)

rule alignment_bam:
  input:
    rules.subset.output.fasta
  output:
    "output/{subset}/{gene}/mapped.bam"
  shell:
    "bealign -r HXB2_{wildcards.gene} {input} {output}"

rule alignment:
  input:
    rules.alignment_bam.output
  output:
    "output/{subset}/{gene}/nuc.fasta"
  shell:
    "bam2msa {input} {output}"

rule protein_alignment:
  input:
    rules.alignment.output
  output:
    "output/{subset}/{gene}/aa.fasta"
  shell:
    "translate {input} > {output}"

rule tn93_distances:
  input:
    rules.alignment.output
  output:
    "output/{subset}/{gene}/tn93.csv"
  shell:
    "tn93 -o {output} -t 1 {input}"

rule kmers:
  input:
    fasta="output/{subset}/{gene}/{char}.fasta"
  output:
    csv="output/{subset}/{gene}/{k}/{char}_kmers.csv"
  run:
    extract_kmer_frequencies(input.fasta, output.csv, wildcards.k)

rule dimensionality_reduction:
  input:
    csv=rules.kmers.output.csv
  output:
    csv="output/{subset}/{gene}/{k}/{char}_{dimension}d_{method}.csv"
  run:
    dimensionality_reduction(input.csv, output.csv, wildcards.dimension, wildcards.method)

rule pairwise_information:
  input:
    tn93=rules.tn93_distances.output[0],
    dr=rules.dimensionality_reduction.output.csv
  output:
    "output/{subset}/{gene}/{k}/pairwise_{char}_{dimension}d_{method}.csv"
  run:
    pairwise_information(input.tn93, input.dr, output[0], wildcards.dimension)

rule plot:
  input:
    script="src/make_plot.py",
    csv=rules.pairwise_information.output
  output:
    png="output/{subset}/{gene}/{k}/{char}_{dimension}d_{method}.png"
  run:
    plot(input.csv[0], output.png)

rule all:
  input:
    expand(
      "output/{subset}/{gene}/{k}/nuc_{dimension}d_{method}.png",
      subset=[125, 250],
      gene=["env", "gag", "pol", "rt"],
      k=[2,3,4,5],
      dimension=[2,3],
      method=["pca", "tsne"]
    )
