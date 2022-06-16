#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Summarise typed isolates')

parser.add_argument('samples', nargs='+', help="List of sample names to summarise", type=str)
parser.add_argument('--species', dest='species', help="Species to summarise", choices=['Escherichia coli', 'Neisseria meningitidis', 'Ecoli', 'Nmen', 'Streptococcus pyogenes', 'Spyo'], type=str, required=True)
parser.add_argument('--qc', dest='qc', help="Path to QC report", type=str, required=True)
parser.add_argument('--fastbaps', dest='fastbaps', help="Path to fastBAPS clusters.csv", type=str, required=True)
parser.add_argument('--mlst', dest="mlst", help="mlst output directory", type=str)
parser.add_argument('--amrfinder', dest="amrfinder", help="AMRfinder output directory", type=str)
parser.add_argument('--vfdb', dest="vfdb", help="ABRicate VFDB output directory", type=str)
parser.add_argument('--emmtyper', dest="emmtyper", help="emmtyper output directory", type=str, default="emmtyper_Spyo")
parser.add_argument('--output', dest='output', help="Output file", type=str, required=True)

args = parser.parse_args()

if (args.species == "Streptococcus pyogenes") or (args.species == "Spyo"):
  args.species = "Spyo"
  assert args.mlst is not None
  assert args.emmtyper is not None

# Store samples as isolate_list
isolate_list = list(args.samples)

# Make output prefix based on timestamp
output_prefix = 'tmp_data/'

list_vfdb_genes = ["fbaA", "fbp54", "grab", "hylP", "ideS/mac", "lmb", "mf2", "mf3", "mf/spd", "scpA", "ska", "slo", "smeZ", "speB", "spec", "ssa"]

def customreadlines(file):
  '''
  Read all lines of a file and close connection.

  Parameters
  ----------
  file : str
    File name to open

  Returns
  -------
  tmp.readlines() : list
    List containing all lines

  '''
  tmp = open(file)
  return tmp.readlines()
  tmp.close()

def parse_table(gene, lines, WT):
  tmp_list = []
  for line in lines:
    if gene in line.split('\t')[5]:
      tmp_list.append(line.split('\t')[5])
  if len(tmp_list) == 0:
    tmp_list.append(WT)
  output_string = '|'.join(tmp_list)
  return output_string

def parse_vfdb(gene, lines):
  tmp_list = []
  for line in lines:
    if gene in line.split('\t')[5]:
      tmp_list.append(line.split('\t')[5])
  if len(tmp_list) == 0:
    output_string = 'None'
  else:
    tmp_list.sort()
    output_string = gene
  for i in tmp_list:
    output_string = output_string + str(i[-1:])
  return output_string

def add_qc_results(full_df, qc_df_path):
  qc_df = pd.read_csv(qc_df_path, sep='\t')
  full_df = qc_df.merge(full_df, on='Isolate', how='right', validate='one_to_one')
  full_df = full_df.fillna('-')
  return full_df

def add_fastbaps(full_df, fastbaps_path):
  fastbaps_df = pd.read_csv(fastbaps_path, names=['index', 'Isolate', 'fastBAPS_lvl1', 'fastBAPS_lvl2'], header=0, dtype={'index': str, 'Isolate': str, 'fastBAPS_lvl1': int, 'fastBAPS_lvl2': int})
  fastbaps_df = fastbaps_df.drop('index', axis=1)
  full_df = full_df.merge(fastbaps_df, on='Isolate', how='left', validate='one_to_one')
  full_df = full_df.fillna(0)
  return full_df

def process_spyo(args):
  full_df = pd.DataFrame()

  for isolate in isolate_list:
    # Get ST
    filepath = output_prefix + args.mlst + '/' + isolate + ".tsv"
    isolate_df = pd.read_csv(filepath, sep='\t')

    # Get emmtyper results
    filepath_emmtyper = output_prefix + args.emmtyper + '/' + isolate + '.tsv'
    emmtyper_df = pd.read_csv(filepath_emmtyper, sep='\t', names=['Isolate_raw', 'nr_blast_hits', 'nr_clusters', 'emm_type', 'emm_gene_position', 'possible_emm_types', 'possible_emm_gene_locations', 'emm_cluster'])
    emmtyper_df['Isolate'] = os.path.splitext(os.path.basename(emmtyper_df['Isolate_raw'][0]))[0]
    emmtyper_df = emmtyper_df.drop(['Isolate_raw', 'nr_blast_hits', 'nr_clusters', 'possible_emm_types', 'possible_emm_gene_locations'], axis=1)

    isolate_df = isolate_df.merge(emmtyper_df, on='Isolate', how='left')

    ### Get VFDB results
    filepath = output_prefix + args.vfdb + '/' + isolate + '.tsv'

    with open(filepath, 'r') as file:
      lines = file.readlines()

    for gene in list_vfdb_genes:
      isolate_df[gene] = parse_table(gene, lines, '-')

    full_df = pd.concat([full_df, isolate_df])

  full_df = add_qc_results(full_df, args.qc)

  full_df = add_fastbaps(full_df, args.fastbaps)

  full_df['fastBAPS_lvl1'] = full_df['fastBAPS_lvl1'].astype(int)
  full_df['fastBAPS_lvl2'] = full_df['fastBAPS_lvl2'].astype(int)

  return full_df


if __name__ == "__main__":
  if args.species == "Spyo":
    full_df = process_spyo(args)

  full_df.to_csv(args.output, index=False)

