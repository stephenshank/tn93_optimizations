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
    "output/{subset}/env.bam"
  shell:
    "bealign -r HXB2_env {input} {output}"

rule alignment:
  input:
    rules.alignment_bam.output
  output:
    "output/{subset}/env_nuc.fasta"
  shell:
    "bam2msa {input} {output}"

rule protein_alignment:
  input:
    rules.alignment.output
  output:
    "output/{subset}/env_aa.fasta"
  shell:
    "translate {input} > {output}"

rule tn93_distances:
  input:
    rules.alignment.output
  output:
    "output/{subset}/env_tn93.csv"
  shell:
    "tn93 -o {output} -t 1 {input}"

rule kmers:
  input:
    fasta="output/{subset}/env_{char}.fasta"
  output:
    csv="output/{subset}/{k}/{char}_kmers.csv"
  run:
    extract_kmer_frequencies(input.fasta, output.csv, wildcards.k)

rule tsne_distances:
  input:
    csv=rules.kmers.output.csv
  output:
    csv="output/{subset}/{k}/{char}_{dimension}d_tsne.csv"
  run:
    tsne_distances(input.csv, output.csv, wildcards.dimension)

rule pairwise_information:
  input:
    rules.tn93_distances.output,
    rules.tsne_distances.output.csv
  output:
    "output/{subset}/{k}/pairwise_{char}_{dimension}d.csv"
  run:
    pairwise_information(wildcards.subset, wildcards.k,  wildcards.char, wildcards.dimension)

rule plot:
  input:
    rules.pairwise_information.output
  output:
    "output/{subset}/{k}/plots/{char}_{dimension}d.png"
  run:
    plot(wildcards.subset, wildcards.k, wildcards.char, wildcards.dimension)

SUBSETS = ["1000"]
NUCLEOTIDE_KS = ["3", "4", "5"]
DIMENSIONS = ["2", "3"]
rule all:
  input:
    expand(
      "output/{subset}/{k}/plots/nuc_{dimension}d.png",
      subset=SUBSETS,
      k=NUCLEOTIDE_KS,
      dimension=DIMENSIONS
    ),
    expand(
      "output/{subset}/2/plots/aa_{dimension}d.png",
      subset=SUBSETS,
      dimension=DIMENSIONS
    )
