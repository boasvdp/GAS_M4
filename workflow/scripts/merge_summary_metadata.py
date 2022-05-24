#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Merge typing summary and metadata.')

parser.add_argument('--metadata', help="Metadata file", type=str, required=True)
parser.add_argument('--summary', help="Summary file", type=str, required=True)
parser.add_argument('--output', help="Output file", type=str, required=True)

args = parser.parse_args()


import pandas as pd
df_summary = pd.read_csv(args.summary)
df_meta = pd.read_csv(args.metadata, sep='\t')
df = df_summary.merge(df_meta, on="Isolate", how="left").fillna('NA')
df.to_csv(args.output, sep='\t', index=False)
