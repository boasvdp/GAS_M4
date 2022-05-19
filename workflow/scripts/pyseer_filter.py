#!/usr/bin/env python3

import argparse

import pandas as pd

def filterdf(df, variants, beta):
  threshold = 0.05 / variants
  filter_col = (df['filter-pvalue'] < threshold)
  if beta == 'both':
    filter_beta = (~df['beta'].isna())
  elif (beta == 'associated') or (beta == 'positive'):
    filter_beta = (df['beta'] > 0)
  elif (beta == 'dissociated') or (beta == 'negative'):
    filter_beta = (df['beta'] < 0)
  df_filtered = df[filter_col & filter_beta]
  return df_filtered

def main(args):
  df = pd.read_csv(args.input, sep='\t', header=0, names=['variant', 'af', 'filter-pvalue', 'lrt-pvalue', 'beta', 'beta-std-err', 'intercept', 'x', 'k-samples', 'nk-samples', 'notes'])
  df = df.fillna('NA')
  print(df.head())
  df = df.sort_values('filter-pvalue', ascending=True)
  df_filtered = filterdf(df, args.variants, args.beta)
  df.to_csv(args.output_sort, sep = '\t',index=False)
  df_filtered.to_csv(args.output_filt, sep = '\t',index=False)

if __name__ == "__main__":
  import argparse
  import pandas as pd
  
  parser = argparse.ArgumentParser(description='Print significant rows from PysEER and sort')

  parser.add_argument('-i', '--input', dest='input', help="Input file", type=str, required=True)
  parser.add_argument('--output-unfiltered', dest='output_sort', help="UNfiltered and sorted output file", type=str, required=True)
  parser.add_argument('--output-filtered', dest='output_filt', help="Filtered and sorted output file", type=str, required=True)
  parser.add_argument('-v', '--variants', dest='variants', help="Number of tested variants", type=int, required=True)
  parser.add_argument('--beta', dest='beta', help='Filter for only positive/negative beta', choices=['positive', 'negative', 'both', 'associated', 'dissociated'], default='both')

  args = parser.parse_args()
  main(args)
