#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description='Convert SNP table to a network graph')

parser.add_argument('-i', '--input', dest='input', help="Input SNP table", type=str, required=True)
parser.add_argument('-o', '--output', dest='output', help="Output Graphml graph", type=str, required=True)
parser.add_argument('-t', '--threshold', dest='threshold', help="SNP calling threshold", type=int, required=True)

args = parser.parse_args()

import pandas as pd
import networkx as nx

df = pd.read_csv(args.input, sep = '\t', header=None, names=["isolate_1", "isolate_2", "SNPs"])

df_filt = df.query('isolate_1 != isolate_2 & SNPs <= @args.threshold')[['isolate_1','isolate_2','SNPs']]
graph = nx.from_pandas_edgelist(df_filt, 'isolate_1', 'isolate_2', 'SNPs')

nx.write_graphml(graph, args.output)
