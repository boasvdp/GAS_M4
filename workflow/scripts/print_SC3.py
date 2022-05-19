#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Print lineage phenotype')

parser.add_argument('--input', help="Input summary file", type=str, required=True)
parser.add_argument('--output', help="Output file", type=str, required=True)

args = parser.parse_args()

df = pd.read_csv(args.input)

df = df[df['QC_decision'] == 'PASS']

df['is_SC3'] = df['M4_type'].map({'SC3': 1}).fillna(0)

df = df[['Isolate', 'is_SC3']]

df['is_SC3'] = df['is_SC3'].astype(int)

df.to_csv(args.output, sep='\t', index=False)
