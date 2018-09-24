#!/usr/bin/env python
# encoding: utf-8
'''

interesting_genes_from_corrs.py


Created by Joan Smith
on 2018-4-2.


Copyright (c) 2018. All rights reserved.

'''

import pandas as pd
import sys
import argparse
import os
import glob

from IPython import embed

sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base

def get_options():
  parser = argparse.ArgumentParser(description='Get cnv directory. csv file')
  parser.add_argument('-i', action='store', dest='corrs_directory')
  parser.add_argument('-f', action='store', dest='interesting_file')
  parser.add_argument('-o', action='store', dest='outdir')
  namespace = parser.parse_args()

  return (namespace.corrs_directory, namespace.interesting_file, namespace.outdir)


def main(argv=None):
  pancan_dir, input_file, outdir = get_options()
  if outdir == None:
    outdir = pancan_dir

  corrs = pd.read_csv(os.path.join(pancan_dir, 'corr_results.csv'))
  interesting_genes = pd.read_csv(input_file, comment='#')
  interesting_genes['Gene'] = '\'' + interesting_genes['Gene']

  corrs['genes'] = '\'' + corrs['Unnamed: 0'].str.split('|', expand=True)[0]
  corrs = corrs.set_index('genes')

  found_corrs = []
  for i, g in interesting_genes.iterrows():
    print g.Gene, g['Cancer Type'], corrs[g['Cancer Type'] + '_' + g.Gene].loc[g.Gene]
    found_corrs.append(corrs[g['Cancer Type'] + '_' + g.Gene].loc[g.Gene])

  interesting_genes['Corrs'] = found_corrs
  interesting_genes.to_csv(os.path.join(outdir, 'interesting_genes.csv'))



if __name__ == "__main__":
  main()
