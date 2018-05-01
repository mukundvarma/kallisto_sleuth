# Toy pipeline for E.Coli analysis using snakemake
# Genome used: E.Coli K12 assembly, can be downloaded from UCSC / Ensembl
# ftp://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/
# Nods to github/slowkow and github/saketrc for snakemake/kallisto tutorials

from os.path import join, basename, dirname
from subprocess import check_output

os.system('mkdir -p reference/kallisto')
os.system('mkdir -p logs')
os.system('mkdir -p outs/objects')
os.system('mkdir -p outs/counts')
os.system('mkdir -p outs/plots')

SAMPLES = ['ecoli_state1_rep1',
  'ecoli_state1_rep2',
  'ecoli_state2_rep1',
  'ecoli_state2_rep2']

rule all:
  input:
    expand('ecoli_cdna/{sample}.fastq.gz', sample = SAMPLES),
    expand('outs/counts/{sample}/abundance.tsv', sample = SAMPLES),
    'metadata.tsv',
    'outs/objects/sleuth.object.RDS'
    
rule kallisto_index:
  input:
    ref = 'reference/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa'
  output:
    index = join('reference', 'kallisto', 'ecoli_kallisto_index.fa')
  log:
    'logs/kallisto.index.log'
  run:
    shell('kallisto index'
      ' --index={output.index}'
      ' --make-unique'
      ' {input.ref}'
      ' >> {log} 2>&1')

rule quantify:
  input:
    fq1 = 'ecoli_cdna/{sample}.fastq.gz',
    fq2 = 'ecoli_cdna/{sample}_2.fastq.gz',
    index = rules.kallisto_index.output.index
  output:
    'outs/counts/{sample}/abundance.tsv'
  params:
    outdir = 'outs/counts/{sample}'
  threads:
    4
  resources:
    mem = 4000
  run:
    shell('kallisto quant --index={input.index}'
        ' --single'
        ' --fragment-length=100 '
        ' --sd=10'
        ' --threads=4'
        ' --seed=13'
#        ' --plaintext'
        ' --output-dir={params.outdir} -b 100 {input.fq1} {input.fq2}')    

rule sleuth:
  input:
    'metadata.tsv'
  output:
    'outs/objects/sleuth.differential_expression.tsv',
    'outs/objects/sleuth.object.RDS'
  run:
     shell('Rscript runSleuth.R')
