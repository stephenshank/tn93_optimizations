from src import *

rule paper:
  input:
    "paper/Draft.tex"
  output:
    "paper/Draft.pdf"
  shell:
    "pdflatex -output-directory=paper {input}"

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
    "bealign -K -r HXB2_{wildcards.gene} {input} {output}"

rule alignment_with_reference:
  input:
    rules.alignment_bam.output
  output:
    "output/{subset}/{gene}/nuc_ref.fasta"
  shell:
    "bam2msa {input} {output}"

rule alignment:
  input:
    rules.alignment_with_reference.output
  output:
    "output/{subset}/{gene}/nuc.fasta"
  run:
    trim_first(input[0], output[0])

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
    script="src/features.py",
    fasta="output/{subset}/{gene}/{char}.fasta"
  output:
    csv="output/{subset}/{gene}/{k}/{char}_kmers.csv"
  run:
    extract_kmer_frequencies(input.fasta, output.csv, wildcards.k)

rule windows:
  input:
    fasta=rules.alignment_with_reference.output[0]
  output:
    csv="output/{subset}/{gene}/window_{length}.csv"
  run:
    window(input.fasta, output.csv, wildcards.length)

rule concat:
  input:
    kmers=rules.kmers.output.csv,
    windows=rules.windows.output.csv
  output:
    csv="output/{subset}/{gene}/{k}/concat_{char}_{length}.csv"
  run:
    concat(input.kmers, input.windows, output.csv)

rule dimensionality_reduction:
  input:
    script="src/run_dimensionality_reduction.py",
    csv=rules.concat.output.csv
  output:
    csv="output/{subset}/{gene}/{k}/{char}_{length}_{dimension}d_{method}.csv"
  run:
    dimensionality_reduction(input.csv, output.csv, wildcards.dimension, wildcards.method)

rule pairwise_information:
  input:
    script="src/make_pairwise_information.py",
    tn93=rules.tn93_distances.output[0],
    dr=rules.dimensionality_reduction.output.csv
  output:
    "output/{subset}/{gene}/{k}/pairwise_{char}_{length}_{dimension}d_{method}.csv"
  run:
    pairwise_information(input.tn93, input.dr, output[0], wildcards.dimension)

rule plot:
  input:
    script="src/output.py",
    csv=rules.pairwise_information.output
  output:
    png="output/{subset}/{gene}/{k}/{char}_{length}_{dimension}d_{method}_{distance}.png"
  run:
    plot(input.csv[0], output.png, wildcards.distance)

rule small:
  input:
    expand(
      "output/{subset}/{gene}/{k}/nuc_{dimension}d_{method}_{distance}.png",
      subset=[125, 250],
      gene=["env", "gag", "pol", "rt"],
      k=[2,3,4,5],
      dimension=[2,3],
      method=["pca", "tsne"],
      distance=['Euclidean', 'L1']
    )
 
TABLE_SUBSETS = range(250, 2001, 250)
rule table:
  input:
    script="src/output.py",
    pairwise=expand(
      "output/{subset}/{{gene}}/{{k}}/pairwise_{{char}}_{{length}}_{{dimension}}d_{{method}}.csv",
      subset=TABLE_SUBSETS
    )
  output:
    "tables/{gene}/{k}/{char}_{length}_{dimension}d_{method}.csv"
  run:
    table(output[0], TABLE_SUBSETS)

