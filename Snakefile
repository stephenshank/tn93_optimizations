rule subset:
  input:
    "input/HIV-LANL-unaligned.fasta"
  output:
    "output/HIV-LANL-unaligned_subset-{subset}.fasta"
  shell:
    "python subset.py -i {input} -o {output} -n {wildcards.subset}"

rule alignment:
  input:
    "output/HIV-LANL-unaligned_subset-{subset}.fasta"
  output:
    "output/HIV-LANL_subset-{subset}.fasta"
  shell:
    "mafft {input} > {output}"

rule tn93_distances:
  input:
    "output/HIV-LANL_subset-{subset}.fasta"
  output:
    "output/HIV-LANL_subset-{subset}_tn93.csv"
  shell:
    "tn93 -o {output} -t 1 {input}"

rule kmers:
  input:
    "output/HIV-LANL-unaligned_subset-{subset}.fasta"
  output:
    "output/HIV-LANL_subset-{subset}_{k}-mers.csv"
  shell:
    "python extract_kmer_frequencies.py -i {input} -o {output} -k {wildcards.k}"

rule tsne_distances:
  input:
    "output/HIV-LANL_subset-{subset}_{k}-mers.csv"
  output:
    "output/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim_tsne.csv"
  shell:
    "python run_tsne.py -i {input} -o {output} -d {wildcards.dimension}"

rule pairwise_information:
  input:
    "output/HIV-LANL_subset-{subset}_tn93.csv",
    "output/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim_tsne.csv"
  output:
    "output/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.csv"
  shell:
    "python make_pairwise_information.py -s {wildcards.subset} -k {wildcards.k} -d {wildcards.dimension}"

rule plot:
  input:
    "output/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.csv"
  output:
    "plots/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.png"
  shell:
    "python make_plot.py -s {wildcards.subset} -k {wildcards.k} -d {wildcards.dimension}"

SUBSETS = ["1000"]
KS = ["3", "4", "5"]
DIMENSIONS = ["2", "3"]
rule all:
  input:
    expand("plots/HIV-LANL_subset-{subset}_{k}-mers_{dimension}-dim.png", subset=SUBSETS, k=KS, dimension=DIMENSIONS)
