#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-07.

Use the interesting genes raw data, and CNAs.csv to count patients with
both mutations and cn variations in the same gene.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np

import glob
import argparse
import sys
import os
import pdb
from multiprocessing import Pool

sys.path.append('../common/')
import mutation_base
import utilities as util
import analysis


CNA_CHANGE_CUTOFF = 0.3
DEEP_CNA_CHANGE_CUTOFF = 1

def get_options():
  parser = argparse.ArgumentParser(description='Find patients with both CNAs and mutations in the same interesting genes')
  parser.add_argument('-i', action='store', dest='interesting_genes_raw')
  parser.add_argument('-f', action='store', dest='interesting_genes_file')
  parser.add_argument('-o', action='store', dest='output_dir', default='.')
  ns = parser.parse_args()

  return (ns.interesting_genes_raw, ns.interesting_genes_file,
          ns.output_dir)

def main():
  interesting_genes_folder, interesting_genes_file, outdir = get_options()
  interesting_genes = pd.read_csv(interesting_genes_file, index_col=[1,2])

  raw_data = os.listdir(interesting_genes_folder)
  raw_data = util.remove_extraneous_files(raw_data)

  results = pd.DataFrame()
  for f in raw_data:
    cancer_type, gene = f.split('_')[0:2]
    change_type = interesting_genes.loc[cancer_type, gene].values[0]
    print cancer_type, gene, change_type

    ig_data = pd.read_csv(os.path.join(interesting_genes_folder, f), index_col=0)

    mutation_id = gene + '_mutations'
    if change_type == 'Deletion':
      ig_data['change'] = ig_data['\'' + gene] <= -CNA_CHANGE_CUTOFF
      ig_data['deep_change'] = ig_data['\'' + gene] <= -DEEP_CNA_CHANGE_CUTOFF
    if change_type == 'Amplification':
      ig_data['change'] = (ig_data['\'' + gene] >= CNA_CHANGE_CUTOFF)
      ig_data['deep_change'] = (ig_data['\'' + gene] >= DEEP_CNA_CHANGE_CUTOFF)


    results = results.append({
                    'Gene': gene,
                    'Cancer Type': cancer_type,
                    'Patinet Count': len(ig_data),
                    'CN Change Count': ig_data['change'].sum(),
                    'CN Deep Change Count': ig_data['deep_change'].sum(),
                    'Mutation Count': ig_data[mutation_id].sum(),
                    'Change + Mutation Count': (ig_data['change'] & ig_data[mutation_id]).sum(),
                    'Deep Change + Mutation Count': (ig_data['deep_change'] & ig_data[mutation_id]).sum()
    }, ignore_index=True)


  results = results.set_index(['Cancer Type', 'Gene'])
  results.to_csv(os.path.join(outdir, 'mutation_and_amplification_overlap.csv'))

if __name__ == "__main__":
  main()
