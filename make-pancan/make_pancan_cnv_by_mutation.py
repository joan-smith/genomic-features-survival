#!/usr/bin/env python
# encoding: utf-8
'''
make pancan_cnv_by_mutation.py

Created by Joan Smith
on 2017-7-9.

Given a set of zscore files for each cancer type, make a pancan file, including stouffers

Copyright (c) 2017 . All rights reserved.
'''

import pandas as pd
import numpy as np
import getopt
import sys
import os

sys.path.append('../common/')
import analysis


CANCER_TYPES = ['BLCA', 'BRCA', 'COADREAD', 'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'LUAD',
                'LUSC', 'OV', 'PAAD', 'PRAD', 'SARC', 'SKCM', 'STES', 'UCEC']

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'hs:i:o:mgc',
    ['help', 'outdir=', 'suffix=', 'indir='])
  except getopt.error, msg:
    usage()

  suffix = None
  indir = None
  outdir = None

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-i', '--indir'):
      indir = value
    if option in ('-s', '--suffix'):
      suffix = value

  if not outdir:
    outdir = indir

  return outdir, indir, suffix

def make_path(indir, cancer, suffix):
  if os.path.isdir(os.path.join(indir, cancer)):
    return os.path.join(indir, cancer, cancer + suffix)
  else:
    return os.path.join(indir, cancer + suffix)

def make_pancan_df(outdir, indir, suffix):
  pancan = pd.DataFrame()
  zscore_header = 'zscore'
  outfile = 'cnv_by_mutation_pancan.csv'

  cancer_types = CANCER_TYPES

  for cancer in cancer_types:
    make_path(indir, cancer, suffix)
    path = make_path(indir, cancer, suffix)
    df = pd.read_csv(path, index_col=0, na_values=[' NA']) # sometimes rpy2 gives back this monstrosity of a NaN value.
    if df.shape[0] > 0:
      if '\'' not in df.index[0]:
        # ADD ' to index
        df = df.reset_index()
        print 'ADD apostrophe'
        df['gene'] = '\'' + df['gene']
        df = df.set_index('gene')
      df = df.rename(columns = lambda x: cancer + ' ' + x.strip())
      pancan = pd.concat((pancan, df), axis=1, verify_integrity=True)

  pancan_df = pd.DataFrame(pancan)

  pancan_df = pancan_df.astype(float)
  print pancan_df.filter(regex=('.* mutated zscore'))
  pancan_df['mutated stouffer unweighted'] = analysis.stouffer_unweighted(pancan_df.filter(regex=('.* mutated zscore')))
  pancan_df['non-mutated stouffer unweighted'] = analysis.stouffer_unweighted(pancan_df.filter(regex=('.* non-mutated zscore')))
  outpath = os.path.join(outdir, outfile)
  pancan_df.to_csv(outpath, index_label='gene')

def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, indir, suffix = get_options(argv)
    make_pancan_df(outdir, indir, suffix)

if __name__ == "__main__":
  main()
