#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Write iTOL colorstrip files for PySEER hits')

parser.add_argument('-i', '--input', dest='input', help="Input file with PySEER hits", type=str, required=True)
parser.add_argument('-o', '--output', dest='outdir', help="Output directory", type=str, required=True)

args = parser.parse_args()

if not os.path.exists(args.outdir):
  os.makedirs(args.outdir)

df = pd.read_csv(args.input, sep='\t')

header_top = ['DATASET_COLORSTRIP', 'SEPARATOR COMMA', 'BORDER_WIDTH,0.5', 'COLOR,#719c4e']
header_middle = ['DATASET_LABEL', 'LEGEND_TITLE']
header_bottom = ['LEGEND_COLORS,#000000,#ffffff', 'LEGEND_LABELS,present,absent', 'LEGEND_SHAPES,1,1', 'MARGIN,5', 'STRIP_WIDTH,50', 'DATA']

for index, row in df.iterrows():
  variant_name = str(row['variant'])
  k_samples = row['k-samples'].split(',')
  nk_samples = row['nk-samples'].split(',')
  outputfile = ''.join(['itol_', variant_name, '.txt'])
  outputpath  = '/'.join([args.outdir, outputfile])
  with open(outputpath, 'w+') as file:
    for line in header_top:
      file.write(line + '\n')
    for line in header_middle:
      writestring = ','.join([line, variant_name])
      file.write(writestring + '\n')
    for line in header_bottom:
      file.write(line + '\n')
    for sample in k_samples:
      writestring = ','.join([str(sample), '#000000', 'present'])
      file.write(writestring + '\n')
    for sample in nk_samples:
      writestring = ','.join([str(sample), '#ffffff', 'absent'])
      file.write(writestring + '\n')
