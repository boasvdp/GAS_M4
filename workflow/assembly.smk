# This defines the file containing configurations
# These values can be accessed through config[entry]
# See https://yaml.org/ for explanation of the YAML format
configfile: "workflow/config/config.yaml"

## Identify sample names from the R1 read files in the input folder
## The script check.py should have already checked whether all samples have a single R1 and R2 file associated with it and whether the input folder does not contain other files
#(SAMPLE_NAMES, SAMPLE_NUMBERS) = glob_wildcards("input/{samplename}_{samplenumber}_L001_R1_001.fastq.gz")
#
## Zip and join together all sample names and numbers. This will be used later on to find files to remove 
#SAMPLE_NAMES_NUMBERS = ['_'.join([samplename, samplenumber]) for samplename, samplenumber in zip(SAMPLE_NAMES, SAMPLE_NUMBERS)]

SAMPLE_NAMES, = glob_wildcards("reads/{samplename}_L001_R1_001.fastq.gz")

# Define the desired output of the pipeline 
rule all:
  input:
#    expand("tmp_data/fastqc_pre_out/{sample}_{read}_fastqc.zip", sample=SAMPLE_NAMES, read = ['R1', 'R2']),
#    expand("tmp_data/fastqc_post_out/{sample}_{read}_fastqc.zip", sample=SAMPLE_NAMES, read = ['R1', 'R2']),
#    expand("tmp_data/kraken_out/{sample}.txt", sample=SAMPLE_NAMES),
#    expand("tmp_data/quast_out/{sample}", sample=SAMPLE_NAMES),
#    expand("backup/{sample}.tar.gz", sample=SAMPLE_NAMES),
    expand("output/summary.csv"),

## Remove illumina sample numbers from read file names and copy to scratch
#rule move_to_scratch:
#  output:
#    fw = "reads/{samplename}_L001_R1_001.fastq.gz",
#    rv = "reads/{samplename}_L001_R2_001.fastq.gz",
#  log:
#    "slurm/snakemake_logs/move_to_scratch/{samplename}.log"
#  shell:
#    """
#    mkdir -p /scratch/reads
#    cp input/{wildcards.samplename}_S*_L001_R1_001.fastq.gz {output.fw}
#    cp input/{wildcards.samplename}_S*_L001_R2_001.fastq.gz {output.rv}
#    """
#
# Run FastQC before read trimming on R1 reads. This uses a Snakemake wrapper
rule fastqc_pre_R1:
  input:
    "reads/{sample}_L001_R1_001.fastq.gz"
  output:
    html = "tmp_data/fastqc_pre_out/{sample}_R1.html",
    zip = "tmp_data/fastqc_pre_out/{sample}_R1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/fastqc_pre/{sample}_R1.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

# Run FastQC before read trimming on R2 reads. This uses a Snakemake wrapper
rule fastqc_pre_R2:
  input:
    "reads/{sample}_L001_R2_001.fastq.gz"
  output:
    html = "tmp_data/fastqc_pre_out/{sample}_R2.html",
    zip = "tmp_data/fastqc_pre_out/{sample}_R2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/fastqc_pre/{sample}_R2.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

# Trim and filter reads using Trimmomatic, executed as a Snakemake wrapper
# See http://www.usadellab.org/cms/?page=trimmomatic for explanation of the Trimmomatic settings
# ILLUMINACLIP:all_paired.fa:2:30:10 Removes Illumina adapters based on FASTA file all_paired.fa. Use 2 seed mismatches, a palindrome clip threshold of 30 and a simple clip threshold of 10
# LEADING:3 Remove low quality bases from the beginning. Minimum quality needed is 3. Special Illumina "low quality segment" regions are marked with quality score 2
# TRAILING:3 Remove low quality bases from the end. Minimum quality needed is 3. Special Illumina "low quality segment" regions are marked with quality score 2
# SLIDINGWINDOW:4:20 Perform sliding window trimming, with a window size of 4 bp and a required quality of 20
# MINLEN:36 Remove reads shorter than 36 bp
rule trimmomatic_pe:
  input:
    r1 = "reads/{sample}_L001_R1_001.fastq.gz",
    r2 = "reads/{sample}_L001_R2_001.fastq.gz"
  output:
    r1 = "output/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    r2 = "output/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
    r1_unpaired = "tmp_data/trimmed/{sample}.1.unpaired.fastq.gz",
    r2_unpaired = "tmp_data/trimmed/{sample}.2.unpaired.fastq.gz",
  log:
    "slurm/snakemake_logs/trimmomatic/{sample}.log"
  params:
    # list of trimmers (see manual)
    trimmer=["ILLUMINACLIP:workflow/db/all_paired.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"],
    # optional parameters
    extra="",
    compression_level="-9"
  threads:
    16
  wrapper:
    "0.78.0/bio/trimmomatic/pe"

# Run FastQC after read trimming on R1 reads. This uses a Snakemake wrapper
rule fastqc_post_R1:
  input:
    "output/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
  output:
    html = "tmp_data/fastqc_post_out/{sample}_R1.html",
    zip = "tmp_data/fastqc_post_out/{sample}_R1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/fastqc_post/{sample}_R1.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

# Run FastQC after read trimming on R2 reads. This uses a Snakemake wrapper
rule fastqc_post_R2:
  input:
    "output/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    html = "tmp_data/fastqc_post_out/{sample}_R2.html",
    zip = "tmp_data/fastqc_post_out/{sample}_R2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/fastqc_post/{sample}_R2.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

# Assemble trimmed data into a draft genome. Only the genome is retained currently.
# The bash script estimates genome size using mash. Shovill tries this as well using kmc but this fails too often in our setup.
# Several parameters are set in the configuration file. Standard settings are to only retain contigs >=500 bp in length,
# subsample reads to estimated 100X depth, use 64 Gb RAM, use SPAdes as assembler and $TMPDIR as tmpdir. 
rule shovill:
  input:
    fw = "output/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "output/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    assembly = "output/genomes/{sample}.fasta",
    shovill = directory("tmp_data/shovill_out/{sample}"),
  conda:
    "envs/shovill.yaml"
  params:
    minlen = config["shovill"]["minlen"],
    ram = config["shovill"]["ram"],
    depth = config["shovill"]["depth"],
    assembler = config["shovill"]["assembler"],
    gsize = config["shovill"]["gsize"],
  log:
    "slurm/snakemake_logs/shovill/{sample}.log"
  threads: 16
  shell:
    """
    shovill --assembler {params.assembler} --outdir {output.shovill} --depth {params.depth} --gsize {params.gsize} --cpus {threads} --ram {params.ram} --minlen {params.minlen} --R1 {input.fw} --R2 {input.rv} 2>&1>{log}
    cp {output.shovill}/contigs.fa {output.assembly}
    """

# Estimate actual read depth. Map trimmed reads back onto assembly to extract depth using bedtools. All data in this step is piped so no intermediate files remain.
# Loosely based on http://thegenomefactory.blogspot.com/2018/10/a-unix-one-liner-to-call-bacterial.html 
rule coverage:
  input:
    fw = "output/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "output/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
    assembly = "output/genomes/{sample}.fasta",
  output:
    "tmp_data/coverage_out/{sample}.txt"
  conda:
    "envs/coverage.yaml"
  params:
    minimap_x = "sr"
  log:
    "slurm/snakemake_logs/coverage/{sample}.log"
  threads: 16
  shell:
    """
    minimap2 -a -x {params.minimap_x} -t {threads} {input.assembly} {input.fw} {input.rv} | samtools sort -l 0 --threads {threads} | bedtools genomecov -d -ibam stdin | awk '{{t += $3}} END {{print t/NR}}' 1>{output} 2>{log}
    """

# Assess assembly quality (without reference genome).
# Version 4.* of Quast is defined in the conda environment file as Quast version 5.* fails to install on some occassions and the most recent version does not offer advantages in our use case
rule quast:
  input:
    "output/genomes/{sample}.fasta"
  output:
    directory("tmp_data/quast_out/{sample}")
  conda:
    "envs/quast.yaml"
  log:
    "slurm/snakemake_logs/quast/{sample}.log"
  threads: 8
  shell:
    """
    quast --threads {threads} -o {output} {input} 2>&1>{log}
    """

# Assess which species are present in the trimmed data.
# Currently this uses a MiniKraken database, but plan to compose a Reflab-specific database with relevant organisms
rule kraken2:
  input:
    fw = "output/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "output/trimmed/{sample}_L001_R2_001_corrected.fastq.gz"
  output:
    report = "tmp_data/kraken_out/{sample}.txt"
  conda:
    "envs/kraken.yaml"
  params:
    general = config["kraken"]["general"],
    db = config["kraken"]["db"]
  log:
    "slurm/snakemake_logs/kraken2/{sample}.log"
  threads: 8
  shell:
    """
    kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
    """

# Identify MLST profile of assembled genome.
rule mlst:
  input:
    "output/genomes/{sample}.fasta"
  output:
    "tmp_data/mlst_out/{sample}.txt"
  log:
    "slurm/snakemake_logs/mlst/{sample}/log"
  conda:
    "envs/mlst.yaml"
  threads: 1
  shell:
    """
    mlst {input} 1>{output} 2>{log}
    """

## Combine all FastQC reports into a multiqc report
#rule multiqc_fastqc:
#  input:
#    expand("tmp_data/fastqc_pre_out/{sample}_{read}_fastqc.zip", sample=SAMPLE_NAMES, read = ['R1', 'R2']),
#    expand("tmp_data/fastqc_post_out/{sample}_{read}_fastqc.zip", sample=SAMPLE_NAMES, read = ['R1', 'R2']),
#  output:
#    "tmp_data/fastqc_report.html"
#  log:
#    "slurm/snakemake_logs/multiqc_fastqc.log"
#  wrapper:
#    "0.78.0/bio/multiqc"
#
## Combine Quast and Kraken reports into a single multiqc report
#rule multiqc_kraken:
#  input:
#    expand("tmp_data/kraken_out/{sample}.txt", sample=SAMPLE_NAMES),
#    expand("tmp_data/quast_out/{sample}", sample=SAMPLE_NAMES)
#  output:
#    "tmp_data/kraken_quast_report.html"
#  log:
#    "slurm/snakemake_logs/multiqc_kraken_quast.log"
#  wrapper:
#    "0.78.0/bio/multiqc"
#
## Copy multiqc reports to the output folder
#rule copy_qc_reports:
#  input:
#    "tmp_data/fastqc_report.html",
#    "tmp_data/kraken_quast_report.html",
#  output:
#    directory("output/qc_reports")
#  threads: 1
#  shell:
#    """
#    mkdir -p {output}
#    cp {input} {output}
#    """

# Find all files to backup and create a tarball with a timestamp
rule backup_data:
  input:
    genome = "output/genomes/{sample}.fasta",
    kraken = "tmp_data/kraken_out/{sample}.txt",
    quast = "tmp_data/quast_out/{sample}",
    fastqc_pre_R1 = "tmp_data/fastqc_pre_out/{sample}_R1.html",
    fastqc_pre_R2 = "tmp_data/fastqc_pre_out/{sample}_R2.html",
    fastqc_post_R1 = "tmp_data/fastqc_post_out/{sample}_R1.html",
    fastqc_post_R2 = "tmp_data/fastqc_post_out/{sample}_R2.html",
    mlst = "tmp_data/mlst_out/{sample}.txt",
    coverage = "tmp_data/coverage_out/{sample}.txt",
  output:
    "backup/{sample}.tar.gz"    
  params:
    sample = "{sample}",
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/quast
    mkdir -p tmp_data/.backup/{params.sample}/fastqc_pre
    mkdir -p tmp_data/.backup/{params.sample}/fastqc_post
    mkdir -p tmp_data/.backup/{params.sample}/kraken
    mkdir -p tmp_data/.backup/{params.sample}/coverage
    mkdir -p tmp_data/.backup/{params.sample}/mlst
    cp {input.fastqc_pre_R1} tmp_data/.backup/{params.sample}/fastqc_pre
    cp {input.fastqc_post_R1} tmp_data/.backup/{params.sample}/fastqc_post
    cp {input.fastqc_pre_R2} tmp_data/.backup/{params.sample}/fastqc_pre
    cp {input.fastqc_post_R2} tmp_data/.backup/{params.sample}/fastqc_post
    cp -r {input.quast} tmp_data/.backup/{params.sample}/quast
    cp {input.kraken} tmp_data/.backup/{params.sample}/kraken
    cp {input.genome} tmp_data/.backup/{params.sample}
    cp {input.mlst} tmp_data/.backup/{params.sample}/mlst
    cp {input.coverage} tmp_data/.backup/{params.sample}/coverage
    tar zcvf backup/{params.sample}.tar.gz -C tmp_data/.backup {params.sample}
    rm -rf tmp_data/.backup/{params.sample}
    """

# Summarise the results in a csv file that will be included in the backup
rule summary:
  input:
    kraken = expand("tmp_data/kraken_out/{sample}.txt", sample=SAMPLE_NAMES),
    quast = expand("tmp_data/quast_out/{sample}", sample=SAMPLE_NAMES),
#    fastqc_pre = expand("tmp_data/fastqc_pre_out/{sample}_{read}_fastqc.zip", sample=SAMPLE_NAMES, read = ['R1', 'R2']),
#    fastqc_post = expand("tmp_data/fastqc_post_out/{sample}_{read}_fastqc.zip", sample=SAMPLE_NAMES, read = ['R1', 'R2']),
    mlst = expand("tmp_data/mlst_out/{sample}.txt", sample=SAMPLE_NAMES),
    coverage = expand("tmp_data/coverage_out/{sample}.txt", sample=SAMPLE_NAMES),
  output:
    "output/summary.csv"
  conda:
    "envs/python.yaml"
  params:
    samples = SAMPLE_NAMES
  log:
    "slurm/snakemake_logs/summary.log"
  threads: 1
  shell:
    """
    python workflow/scripts/summary.py {params.samples} 1>{output} 2>{log}
    """
