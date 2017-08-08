#!/usr/bin/env python
# encoding: utf-8
'''
fdr_correction.py

Created by Joan Smith
on 2017-8-4.

Calculate fdr and corrected pvalues for zscore calculations.
'''

import pprint
import argparse
import sys
import os
import fnmatch
import pandas as pd
import statsmodels.stats.multitest as smm
from collections import defaultdict
import scipy.stats

sys.path.append('../common/')
import utilities as util

def get_options():
  parser = argparse.ArgumentParser(description='Calculate fdr and corrected pvalues for zscores')
  parser.add_argument('-i', action='store', dest='input')
  parser.add_argument('-p', action='store', dest='pattern')
  parser.add_argument('-n', action='store', dest='outname', default='pancan_corrected_pvalues.csv')
  parser.add_argument('-e', action='store', dest='exclude', default='')
  parser.add_argument('-s', action='store_true', dest='subdirectories')
  namespace = parser.parse_args()
  return (namespace.input, namespace.pattern, namespace.outname,
          namespace.exclude, namespace.subdirectories)


def collect_matches(indir, pattern, exclude, subdirectories):
  matches = defaultdict(list)
  for root, dirnames, filenames in os.walk(indir):
    for filename in fnmatch.filter(filenames, pattern):
      if len(exclude) == 0 or not exclude in filename:
        match_name = os.path.join(root, filename)
        directory = os.path.dirname(match_name)
        if subdirectories:
          directory = os.path.dirname(directory)
        matches[directory].append(match_name)
  return matches

def stouffer_fdr(pancan_f):
  print pancan_f
  df = pd.read_csv(pancan_f, index_col=0)
  df = df.astype(float)
  df = df.dropna(subset=['stouffer unweighted'])
  df['pvalue'] = scipy.stats.norm.sf(abs(df['stouffer unweighted']))*2
  rejected, corrected, _, _ = smm.multipletests(df['pvalue'], method='fdr_bh')
  outdf = df[['stouffer unweighted', 'pvalue']].copy()
  outdf['rejected?'] = rejected
  outdf['corrected'] = corrected
  return outdf

def single_zscore_file_fdr(f):
  df = pd.read_csv(f, index_col=0)
  df = df.astype(float)
  df = df.dropna(subset=['pvalue'])
  rejected, corrected, _, _ = smm.multipletests(df.pvalue, method='fdr_bh')
  outdf = df[['pvalue', 'zscore']].copy()
  outdf['rejected?'] = rejected
  outdf['corrected'] = corrected
  if not '\'' in outdf.index.values[0]:
    outdf.index = '\'' + outdf.index.values
  return outdf

def pancan_fdr(directory, files, outname):
  pancan_fdr = pd.DataFrame()

  for f in files:
    cancer_type = util.get_cancer_type(f).split('_')[0]
    print cancer_type
    cancer_type_df = single_zscore_file_fdr(f)
    cancer_type_df = cancer_type_df.add_prefix(cancer_type + ' ')
    pancan_fdr = pd.concat((pancan_fdr, cancer_type_df), axis=1)

  stouffer_fdr_df = stouffer_fdr(os.path.join(directory, 'pancan.csv'))
  stouffer_fdr_df = stouffer_fdr_df.add_prefix('pancan ')
  pancan_fdr = pd.concat((pancan_fdr, stouffer_fdr_df), axis=1, verify_integrity=True)

  pancan_fdr.to_csv(os.path.join(directory, outname))

def main():
  pp = pprint.PrettyPrinter(indent=2)
  indir, pattern, outname, exclude, subdirectories = get_options()
  matches = collect_matches(indir, pattern, exclude, subdirectories)

  pp.pprint(matches.keys())
  for directory in matches.keys():
    pancan_fdr(directory, matches[directory], outname)




if __name__ == "__main__":
  main()
