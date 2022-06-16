import sys

configfile: "workflow/config/config.yaml"

def get_all_samples_from_qc_report(qc_report_file):
  tmp_list = []
  full_path_qc_file = ''.join(['backup/', qc_report_file])
  with open(full_path_qc_file) as qc_file:
    lines = qc_file.readlines()[1:]
  for line in lines:
    tmp_list.append(line.split('\t')[0])
  return tmp_list

def get_qc_pass_samples_from_qc_report(qc_report_file):
  tmp_list = []
  full_path_qc_file = ''.join(['backup/', qc_report_file])
  with open(full_path_qc_file) as qc_file:
    lines = qc_file.readlines()
  for line in lines:
    if line.split('\t')[-1].rstrip('\n') == 'PASS':
      tmp_list.append(line.split('\t')[0])
  return tmp_list

SPYO_SAMPLES = get_all_samples_from_qc_report('qc_report_spyo.tsv')
SPYO_SAMPLES_PASS = get_qc_pass_samples_from_qc_report('qc_report_spyo.tsv')

for isolate in ['SRR18923745', 'SRR18923745', 'SRR18923745', 'SRR18923745']:
  if isolate in SPYO_SAMPLES:
    SPYO_SAMPLES.remove(isolate)
  if isolate in SPYO_SAMPLES_PASS:
    SPYO_SAMPLES_PASS.remove(isolate)

for isolate in ['SRR18932517', 'SRR18917485']:
  if isolate in SPYO_SAMPLES_PASS:
    SPYO_SAMPLES_PASS.remove(isolate)

all_output = []

if len(SPYO_SAMPLES) > 0:
  all_output.append("output/Spyo_typing_summary.csv")
  all_output.append(expand("backup/Spyo_typing/{sample}_typing.tar.gz", sample=SPYO_SAMPLES))
  all_output.append("tmp_data/clonalframeml_out")
  all_output.append("tmp_data/fastBAPS_multi_out/clusters.csv")
  all_output.append("tmp_data/panaroo_out")
  all_output.append("tmp_data/snp_dists_out/snps.tsv")
  all_output.append(expand("tmp_data/abricate_emm4/{sample}.tsv", sample=SPYO_SAMPLES))
  all_output.append(expand("tmp_data/abricate_capsule/{sample}.tsv", sample=SPYO_SAMPLES))
  all_output.append(expand("tmp_data/abricate_vfdb/{sample}.tsv", sample=SPYO_SAMPLES))
  all_output.append("tmp_data/iqtree_masked_out")
#  all_output.append(expand("tmp_data/abricate_GAS_M4_vf/{sample}.tsv", sample=SPYO_SAMPLES))
#  all_output.append("tmp_data/pyseer_phenotype.tsv")
#  all_output.append("tmp_data/pyseer_pres/filtered.txt")
#  all_output.append("tmp_data/pyseer_vcf/filtered.txt")
#  all_output.append("tmp_data/pyseer_struct/filtered.txt")
#  all_output.append("tmp_data/itol_pyseer_hits/pres")
  all_output.append("tmp_data/itol_summary_out")
  all_output.append(expand("tmp_data/bakta_extracted_genes_nga_ifs_slo/{sample}", sample=SPYO_SAMPLES_PASS))
  all_output.append("tmp_data/SNP_network/SNP_network.graphml")
  all_output.append(expand("tmp_data/phigaro/{sample}/{sample}.html", sample=SPYO_SAMPLES))

if len(all_output) == 0:
  print("No samples to process", file=sys.stderr)

rule all:
  input:
    all_output

# Streptococcus pyogenes typing
rule MLST_spyo:
  input:
    genome = "output/genomes/{sample}.fasta"
  output:
    json = "tmp_data/MLST_Spyo/{sample}.json"
  threads: 8
  params:
    scheme = config['pubmlst']['spyo_mlst_scheme']
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/MLST_Spyo_PubMLST_API/{sample}.log"
  shell:
    """
    python workflow/scripts/type_pubmlst_api.py --genome {input.genome} --api-url {params.scheme} 1> {output.json} 2>{log}
    """

rule convert_MLST_spyo:
  input:
    "tmp_data/MLST_Spyo/{sample}.json"
  output:
    "tmp_data/MLST_Spyo/{sample}.tsv"
  threads: 1
  conda:
    "envs/python.yaml"
  params:
    loci = config['pubmlst']['spyo_mlst_loci']
  log:
    "slurm/snakemake_logs/MLST_Spyo_convert/{sample}.log"
  shell:
    """
    python workflow/scripts/convert_MLST.py --input {input} --output {output} --loci {params.loci}
    """

rule emmtyper_spyo:
  input:
    genome = "output/genomes/{sample}.fasta"
  output:
    tsv = "tmp_data/emmtyper_Spyo/{sample}.tsv"
  threads: 8
  conda:
    "envs/emmtyper.yaml"
  log:
    "slurm/snakemake_logs/emmtyper_Spyo/{sample}.log"
  shell:
    """
    emmtyper --output-format verbose {input.genome} > {output.tsv} 2>{log}
    """

rule bakta:
  input:
    genome = "output/genomes/{sample}.fasta",
  output:
    gff = "output/bakta_out/{sample}/{sample}.gff3"
  threads: 8
  conda:
    "envs/bakta.yaml"
  log:
    "slurm/snakemake_logs/bakta/{sample}.log"
  params:
    outdir = "output/bakta_out/{sample}",
    prefix = "{sample}",
    locustag = "{sample}",
    db = config["bakta"]["db"],
    min_length = config["bakta"]["min_length"],
  shell:
    """
    bakta --db {params.db} --keep-contig-headers --min-contig-length {params.min_length} --prefix {params.prefix} --output {params.outdir} --genus Streptococcus --species pyogenes --gram '?' --translation-table 11 --locus-tag {params.locustag} --verbose --threads {threads} {input.genome} 2>&1>{log}
    """

rule gunzip_ref:
  input:
    "workflow/references/MGAS5005.gbk.gz"
  output:
    "tmp_data/references/MGAS5005.gbk"
  threads: 16
  shell:
    """
    gunzip -c {input} > {output}
    """

rule snippy_SC3_Spyo:
  input:
    r1 = "output/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    r2 = "output/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
    ref = "workflow/references/HKU360.gbk"
  output:
    directory("tmp_data/snippy_SC3_Spyo/{sample}")
  log:
    "slurm/snakemake_logs/snippy_SC3_Spyo/{sample}.log"
  params:
    general = config["snippy"]["general"]
  threads: 16
  conda: "envs/snippy.yaml"
  shell:
    """
    snippy {params.general} --cpus {threads} --outdir {output} --ref {input.ref} --R1 {input.r1} --R2 {input.r2} 2>&1>{log}
    """

rule compare_SNPs:
  input:
    snippy = "tmp_data/snippy_SC3_Spyo/{sample}",
    ref = "workflow/references/Duke.vcf"
  output:
    "tmp_data/SC3_SNPs_Spyo/{sample}.tsv"
  conda:
    "envs/compare_SNPs.yaml"
  threads: 1
  log:
    "slurm/snakemake_logs/compare_SNPs/{sample}.log"
  shell:
    """
    bgzip --keep --force {input.snippy}/snps.vcf 2>&1>{log}
    bcftools index {input.snippy}/snps.vcf.gz 2>&1>>{log}
    bcftools view -O v -R {input.ref} {input.snippy}/snps.vcf.gz 2>>{log} | vt decompose_blocksub -o - - 2>>{log} | python workflow/scripts/compare_SNPs.py --ref {input.ref} --extended > {output} 2>>{log}
    """

rule AMRFinder:
  input:
    genome = "output/genomes/{sample}.fasta"
  output: 
    "tmp_data/AMRfinder/{sample}.tsv"
  conda:
    "envs/amrfinder.yaml"
  params: 
    organism = "Streptococcus_pyogenes"
  threads: 16
  log:
    "slurm/snakemake_logs/AMRfinder/{sample}.log"
  shell:
    """
    amrfinder -u
    amrfinder --threads 4 --nucleotide {input} --organism {params.organism} --output {output} 2>&1>{log}    
    """

rule phigaro:
  input:
    genome = "output/genomes/{sample}.fasta",
    config = "workflow/config/phigaro.config"
  output:
    html = "tmp_data/phigaro/{sample}/{sample}.html"
  conda:
    "envs/phigaro.yaml"
  params:
    outname = "output/phigaro/{sample}/{sample}",
    general = "--not-open"
  log:
    "slurm/snakemake_logs/phigaro/{sample}.log"
  threads: 15
  shell:
    """
    phigaro --config {input.config} --fasta-file {input.genome} --output {params.outname} --extension html,tsv --print-vogs --threads {threads} --save-fasta 2>&1>{log}
    """

# Find all files to backup and create a tarball with a timestamp
rule backup_data_spyo:
  input:
    SC3_tsv = "tmp_data/SC3_SNPs_Spyo/{sample}.tsv",
    emmtyper_tsv = "tmp_data/emmtyper_Spyo/{sample}.tsv",
    MLST_tsv = "tmp_data/MLST_Spyo/{sample}.tsv",
    MLST_json = "tmp_data/MLST_Spyo/{sample}.json",
    AMRfinder = "tmp_data/AMRfinder/{sample}.tsv",
  output:
    "backup/Spyo_typing/{sample}_typing.tar.gz"
  params:
    sample = "{sample}",
  threads: 16
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/MLST_PubMLST
    mkdir -p tmp_data/.backup/{params.sample}/emmtyper
    mkdir -p tmp_data/.backup/{params.sample}/SC3_SNPs
    mkdir -p tmp_data/.backup/{params.sample}/AMRfinder
    cp {input.emmtyper_tsv} tmp_data/.backup/{params.sample}/emmtyper
    cp {input.MLST_tsv} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.MLST_json} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.SC3_tsv} tmp_data/.backup/{params.sample}/SC3_SNPs
    cp {input.AMRfinder} tmp_data/.backup/{params.sample}/AMRfinder
    tar zcvf {output} -C tmp_data/.backup {params.sample}
    rm -rf tmp_data/.backup
    """

rule typing_summary:
  input:
    expand("tmp_data/SC3_SNPs_Spyo/{sample}.tsv", sample=SPYO_SAMPLES),
    expand("tmp_data/emmtyper_Spyo/{sample}.tsv", sample=SPYO_SAMPLES),
    expand("tmp_data/MLST_Spyo/{sample}.tsv", sample=SPYO_SAMPLES),
    expand("tmp_data/abricate_vfdb/{sample}.tsv", sample=SPYO_SAMPLES),
    qc_report = "backup/qc_report_spyo.tsv",
    fastbaps = "tmp_data/fastBAPS_multi_out/clusters.csv",
  output:
    "output/Spyo_typing_summary.csv"
  conda:
    "envs/python.yaml"
  params:
    samples = SPYO_SAMPLES
  threads: 16
  log:
    "slurm/snakemake_logs/Spyo_typing_summary.log"
  shell:
    """
    python workflow/scripts/typing_summary.py --species Spyo --fastbaps {input.fastbaps} --qc {input.qc_report} --output {output} --mlst MLST_Spyo --emmtyper emmtyper_Spyo --vfdb abricate_vfdb {params.samples} 2>&1>{log}
    """

rule panaroo:
  input:
    expand("output/bakta_out/{sample}/{sample}.gff3", sample=SPYO_SAMPLES_PASS),
  output:
    directory("tmp_data/panaroo_out")
  conda:
    "envs/panaroo.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/panaroo.log"
  shell:
    """
    mkdir -p {output}
    panaroo -i {input} -o {output} --threads {threads} --clean-mode strict -a pan 2>&1>{log}
    """

rule snippy_core:
  input:
    snippy_data = expand("tmp_data/snippy_SC3_Spyo/{sample}", sample=SPYO_SAMPLES_PASS),
    ref = "workflow/references/HKU360.gbk"
  output:
    full = "tmp_data/snippy_core_out/core.full.aln",
    snps = "tmp_data/snippy_core_out/core.aln",
    vcf = "tmp_data/snippy_core_out/core.vcf",
  conda:
    "envs/snippy.yaml"
  params:
    outdir = "tmp_data/snippy_core_out"
  log:
    "slurm/snakemake_logs/snippy_core.log"
  threads: 16
  shell:
    """
    mkdir -p {params.outdir}
    snippy-core --ref {input.ref} {input.snippy_data} 2>&1>{log}
    mv core.aln core.full.aln core.tab core.vcf core.txt core.ref.fa {params.outdir}
    """

rule fastbaps:
  input:
    "tmp_data/maskrc_out/masked_snps.aln",
  output:
    "tmp_data/fastBAPS_multi_out/clusters.csv",
  conda:
    "envs/fastbaps.yaml"
  log:
    "slurm/snakemake_logs/fastbaps.log"
  threads: 16
  shell:
    """
    Rscript workflow/scripts/fastBAPS_multi.R {input} {output} 2>&1>{log}
    """

rule iqtree:
  input:
    core = "tmp_data/snippy_core_out/core.aln",
    fullcore = "tmp_data/snippy_core_out/core.full.aln"
  output:
    directory("tmp_data/iqtree_out")
  conda:
    "envs/iqtree_snp-sites.yaml"
  params:
    prefix = "Spyo"
  log:
    "slurm/snakemake_logs/iqtree.log"
  threads: 16
  shell:
    """
    mkdir -p {output} && cd {output}
    iqtree -fconst $(snp-sites -C ../../{input.fullcore}) -nt AUTO -pre {params.prefix} -m HKY+F -bb 1000 -s ../../{input.core} 2>&1>../../{log}
    if [ -f {params.prefix}.treefile ]; then echo "{output}/{params.prefix}.treefile exists"; else exit 1; fi
    """

rule clonalframeml:
  input:
    tree = "tmp_data/iqtree_out",
    aln = "tmp_data/snippy_core_out/core.full.aln"
  output:
    directory("tmp_data/clonalframeml_out")
  conda:
    "envs/clonalframeml.yaml"
  params:
    prefix = "Spyo_cfml",
    iqtreeprefix = "Spyo"
  log:
    "slurm/snakemake_logs/clonalframeml.log"
  threads: 16
  shell:
    """
    mkdir -p {output} && cd {output}
    ClonalFrameML ../../{input.tree}/{params.iqtreeprefix}.treefile ../../{input.aln} {params.prefix} -show_progress true 2>&1>../../{log}
    if [ -f {params.prefix}.labelled_tree.newick ]; then echo "CFML output exists"; else exit 1; fi
    """

rule maskrc:
  input:
    aln = "tmp_data/snippy_core_out/core.full.aln",
    cfml = "tmp_data/clonalframeml_out"
  output:
    "tmp_data/maskrc_out/masked.aln"
  conda:
    "envs/maskrc.yaml"
  params:
    prefix = "Spyo_cfml"
  log:
    "slurm/snakemake_logs/maskrc.log"
  shell:
    """
    bash workflow/scripts/download_maskrc.sh
    cd {input.cfml}
    python3 ../../workflow/scripts/maskrc-svg.py --aln ../../{input.aln} --out ../../{output} {params.prefix} 2>&1>../../{log}
    """

rule snp_sites:
  input:
    "tmp_data/maskrc_out/masked.aln"
  output:
    "tmp_data/maskrc_out/masked_snps.aln"
  conda:
    "envs/iqtree_snp-sites.yaml"
  log:
    "slurm/snakemake_logs/snp_sites.log"
  threads: 16
  shell:
    """
    snp-sites {input} > {output} 2>{log}
    """

rule iqtree_masked:
  input:
    masked = "tmp_data/maskrc_out/masked.aln"
  output:
    directory("tmp_data/iqtree_masked_out")
  conda:
    "envs/iqtree_snp-sites.yaml"
  params:
    prefix = "Spyo_masked"
  log:
    "slurm/snakemake_logs/iqtree_masked.log"
  threads: 16
  shell:
    """
    mkdir -p {output} && cd {output}
    iqtree -nt AUTO -pre {params.prefix} -m HKY+F -bb 1000 -s ../../{input.masked} 2>&1>../../{log}
    if [ -f {params.prefix}.treefile ]; then echo "{output}/{params.prefix}.treefile exists"; else exit 1; fi
    """

rule snp_dists:
  input:
    fullcore = "tmp_data/snippy_core_out/core.full.aln",
    masked = "tmp_data/maskrc_out/masked.aln"
  output:
    fullcore = "tmp_data/snp_dists_out/snps.tsv",
    masked = "tmp_data/snp_dists_out/snps_masked.tsv"
  conda:
    "envs/snp_dists.yaml"
  log:
    "slurm/snakemake_logs/snp_dists.log"
  threads: 16
  shell:
    """
    snp-dists -m {input.fullcore} > {output.fullcore} 2>{log}
    snp-dists -m {input.masked} > {output.masked} 2>>{log}
    """

rule abricate_capsule:
  input:
    genome = "output/genomes/{sample}.fasta",
    db = "workflow/db/Spyo_capsule"
  output:
    "tmp_data/abricate_capsule/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    datadir = config["abricate"]["datadir"],
    db = config["abricate"]["capsule_db"],
    mincov = config["abricate"]["mincov"],
    minid = config["abricate"]["minid"],
  threads: 2
  log:
    "slurm/snakemake_logs/abricate_capsule/{sample}.log"
  shell:
    """
    abricate --datadir {params.datadir} --db {params.db} --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output} 2>{log}
    """

rule abricate_emm4:
  input:
    genome = "output/genomes/{sample}.fasta",
    db = "workflow/db/emm4"
  output:
    "tmp_data/abricate_emm4/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    datadir = config["abricate"]["datadir"],
    db = config["abricate"]["emm4_db"],
    mincov = config["abricate"]["mincov"],
    minid = config["abricate"]["minid"],
  threads: 2
  log:
    "slurm/snakemake_logs/abricate_emm4/{sample}.log"
  shell:
    """
    abricate --datadir {params.datadir} --db {params.db} --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output} 2>{log}
    """

rule abricate_vfdb:
  input:
    genome = "output/genomes/{sample}.fasta",
  output:
    "tmp_data/abricate_vfdb/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    db = config["abricate"]["vfdb_db"],
    mincov = config["abricate"]["mincov"],
    minid = config["abricate"]["minid"],
  threads: 2
  log:
    "slurm/snakemake_logs/abricate_vfdb/{sample}.log"
  shell:
    """
    abricate --db {params.db} --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output} 2>{log}
    """

rule abricate_GAS_M4_vf:
  input:
    genome = "output/genomes/{sample}.fasta",
    db = "workflow/db/GAS_M4_vf"
  output:
    "tmp_data/abricate_GAS_M4_vf/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    datadir = config["abricate"]["datadir"],
    db = config["abricate"]["GAS_M4_vf_db"],
    mincov = config["abricate"]["mincov"],
    minid = config["abricate"]["minid"],
  threads: 2
  log:
    "slurm/snakemake_logs/abricate_GAS_M4_vf/{sample}.log"
  shell:
    """
    abricate --datadir {params.datadir} --db {params.db} --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output} 2>{log}
    """

rule abricate_GAS_M4_nga_ifs_slo:
  input:
    genome = "output/genomes/{sample}.fasta",
    db = "workflow/db/GAS_M4_nga_ifs_slo"
  output:
    "tmp_data/abricate_GAS_M4_nga_ifs_slo/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    datadir = config["abricate"]["datadir"],
    db = config["abricate"]["GAS_M4_nga_ifs_slo_db"],
    mincov = config["abricate"]["mincov"],
    minid = config["abricate"]["minid"],
  threads: 2
  log:
    "slurm/snakemake_logs/abricate_GAS_M4_vf/{sample}.log"
  shell:
    """
    abricate --datadir {params.datadir} --db {params.db} --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output} 2>{log}
    """

rule convert_phenotype:
  input:
    "output/Spyo_typing_summary.csv"
  output:
    "tmp_data/pyseer_phenotype.tsv"
  conda:
    "envs/python.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/convert_phenotype.log"
  shell:
    """
    python workflow/scripts/print_SC3.py --input {input} --output {output}
    """

rule pyseer_pres:
  input:
    panaroo = "tmp_data/panaroo_out",
    phenotypes = "workflow/private/pyseer_new_lineage.tsv",
  output:
    selected = "tmp_data/pyseer_pres/selected.txt",
    patterns = "tmp_data/pyseer_pres/patterns.txt",
  conda:
    "envs/pyseer.yaml"
  params:
    general = config["pyseer"]["general"]
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_pres.log"
  shell:
    """
    bash workflow/scripts/download_count_patterns.sh
    pyseer {params.general} --phenotypes {input.phenotypes} --pres {input.panaroo}/gene_presence_absence.Rtab --output-patterns {output.patterns} 1> {output.selected} 2>{log}
    """

rule pyseer_vcf:
  input:
    vcf = "tmp_data/snippy_core_out/core.vcf",
    phenotypes = "workflow/private/pyseer_new_lineage.tsv",
  output:
    selected = "tmp_data/pyseer_vcf/selected.txt",
    patterns = "tmp_data/pyseer_vcf/patterns.txt",
  conda:
    "envs/pyseer.yaml"
  params:
    general = config["pyseer"]["general"]
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_vcf.log"
  shell:
    """
    bash workflow/scripts/download_count_patterns.sh
    pyseer {params.general} --phenotypes {input.phenotypes} --vcf {input.vcf} --output-patterns {output.patterns} 1> {output.selected} 2>{log}
    """

rule pyseer_struct:
  input:
    panaroo = "tmp_data/panaroo_out",
    phenotypes = "workflow/private/pyseer_new_lineage.tsv",
  output:
    selected = "tmp_data/pyseer_struct/selected.txt",
    patterns = "tmp_data/pyseer_struct/patterns.txt",
  conda:
    "envs/pyseer.yaml"
  params:
    general = config["pyseer"]["general"]
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_struct.log"
  shell:
    """
    bash workflow/scripts/download_count_patterns.sh
    pyseer {params.general} --phenotypes {input.phenotypes} --pres {input.panaroo}/struct_presence_absence.Rtab --output-patterns {output.patterns} 1> {output.selected} 2>{log}
    """

rule pyseer_filter_pres:
  input:
    selected = "tmp_data/pyseer_pres/selected.txt",
    patterns = "tmp_data/pyseer_pres/patterns.txt",
  output:
    filtered = "tmp_data/pyseer_pres/filtered.txt",
    unfiltered = "tmp_data/pyseer_pres/unfiltered.txt",
    counted_patterns = "tmp_data/pyseer_pres/counted_patterns.txt",
  conda:
    "envs/pyseer.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_filter_pres.log"
  shell:
    """
    VARIANTS=$(python workflow/scripts/count_patterns.py {input.patterns} | tee {output.counted_patterns} | awk 'NR == 1 {{print $2}}')
    python workflow/scripts/pyseer_filter.py -i {input.selected} --output-filtered {output.filtered} --output-unfiltered {output.unfiltered} --variants $VARIANTS
    """

rule pyseer_filter_vcf:
  input:
    selected = "tmp_data/pyseer_vcf/selected.txt",
    patterns = "tmp_data/pyseer_vcf/patterns.txt",
  output:
    filtered = "tmp_data/pyseer_vcf/filtered.txt",
    unfiltered = "tmp_data/pyseer_vcf/unfiltered.txt",
    counted_patterns = "tmp_data/pyseer_vcf/counted_patterns.txt",
  conda:
    "envs/pyseer.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_filter_vcf.log"
  shell:
    """
    VARIANTS=$(python workflow/scripts/count_patterns.py {input.patterns} | tee {output.counted_patterns} | awk 'NR == 1 {{print $2}}')
    python workflow/scripts/pyseer_filter.py -i {input.selected} --output-filtered {output.filtered} --output-unfiltered {output.unfiltered} --variants $VARIANTS
    """

rule pyseer_filter_struct:
  input:
    selected = "tmp_data/pyseer_struct/selected.txt",
    patterns = "tmp_data/pyseer_struct/patterns.txt",
  output:
    filtered = "tmp_data/pyseer_struct/filtered.txt",
    unfiltered = "tmp_data/pyseer_struct/unfiltered.txt",
    counted_patterns = "tmp_data/pyseer_struct/counted_patterns.txt",
  conda:
    "envs/pyseer.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_filter_struct.log"
  shell:
    """
    VARIANTS=$(python workflow/scripts/count_patterns.py {input.patterns} | tee {output.counted_patterns} | awk 'NR == 1 {{print $2}}')
    python workflow/scripts/pyseer_filter.py -i {input.selected} --output-filtered {output.filtered} --output-unfiltered {output.unfiltered} --variants $VARIANTS
    """

rule pyseer_to_itol:
  input:
    pres = "tmp_data/pyseer_pres/filtered.txt",
    vcf = "tmp_data/pyseer_vcf/filtered.txt",
    struct = "tmp_data/pyseer_struct/filtered.txt",
  output:
    pres = directory("tmp_data/itol_pyseer_hits/pres"),
    vcf = directory("tmp_data/itol_pyseer_hits/vcf"),
    struct = directory("tmp_data/itol_pyseer_hits/struct"),
  conda:
    "envs/python.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/pyseer_to_itol.log"
  shell:
    """
    python workflow/scripts/pyseer_hits_itol.py -i {input.pres} -o {output.pres} 2>&1>{log}
    python workflow/scripts/pyseer_hits_itol.py -i {input.vcf} -o {output.vcf} 2>&1>>{log}
    python workflow/scripts/pyseer_hits_itol.py -i {input.struct} -o {output.struct} 2>&1>>{log}
    """

rule extract_genes_abricate_nga_ifs_slo:
  input:
    expand("output/genomes/{sample}.fasta", sample=SPYO_SAMPLES),
    abricate = "tmp_data/abricate_GAS_M4_nga_ifs_slo/{sample}.tsv",
  output:
    directory("tmp_data/extracted_genes_nga_ifs_slo/{sample}"),
  conda:
    "envs/biopython.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/extract_genes_abricate_nga_ifs_slo/{sample}.log"
  shell:
    """
    python workflow/scripts/extract_genes_abricate.py -a {input.abricate} -g output/genomes -o {output} --genecluster --flanking --flanking-bp 2000 2>&1>{log}
    """

rule bakta_extracted_gened_nga_ifs_slo:
  input:
    "tmp_data/extracted_genes_nga_ifs_slo/{sample}"
  output:
    directory("tmp_data/bakta_extracted_genes_nga_ifs_slo/{sample}")
  conda:
    "envs/bakta.yaml"
  threads: 8
  log:
    "slurm/snakemake_logs/bakta_extracted_gened_nga_ifs_slo/{sample}.log"
  params:
    prefix = "{sample}",
    locustag = "{sample}",
    db = config["bakta"]["db"],
  shell:
    """
    bakta --db {params.db} --keep-contig-headers --prefix {params.prefix} --output {output} --genus Streptococcus --species pyogenes --gram '?' --translation-table 11 --locus-tag {params.locustag} --verbose --threads {threads} {input}/{params.prefix}*.out 2>&1>{log}
    """

rule snp_network:
  input:
    normal = "tmp_data/snippy_core_out/core.full.aln",
    masked = "tmp_data/maskrc_out/masked.aln"    
  output:
    normal = "tmp_data/SNP_network/SNP_network.graphml",
    masked = "tmp_data/SNP_network/SNP_network_masked.graphml"
  conda:
    "envs/snp_network.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/snp_network.log"
  params:
  shell:
    """
    python workflow/scripts/produce_SNP_network.py -i {input.normal} -o {output.normal} -t 10 2>&1>{log}
    python workflow/scripts/produce_SNP_network.py -i {input.masked} -o {output.masked} -t 10 2>&1>>{log}
    """

#rule merge_summary_metadata:
#  input:
#    metadata = "workflow/private/metadata_clean.txt",
#    summary = "output/Spyo_typing_summary.csv",
#  output:
#    merged = "workflow/private/summary_metadata.tsv"
#  conda:
#    "envs/python.yaml"
#  log:
#    "slurm/snakemake_logs/merge_summary_metadata.log"
#  threads: 16
#  shell:
#    """
#    python workflow/scripts/merge_summary_metadata.py --metadata {input.metadata} --summary {input.summary} --output {output.merged} 2>&1>{log} 
#    """

rule itolparser:
  input:
    "output/Spyo_typing_summary.csv"
  output:
    directory("tmp_data/itol_summary_out")
  conda:
    "envs/itolparser.yaml"
  threads: 16
  log:
    "slurm/snakemake_logs/itolparser.log"
  shell:
    """
    ../sentinel/sentinel_analysis/itolparser/itolparser -i {input} -o {output} --ignore emm_gene_position,variant_match --continuous exact_match,non_match -d ','
    """
