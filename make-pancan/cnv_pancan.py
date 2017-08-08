#!/usr/bin/env python
# encoding: utf-8
'''
make pancan.py

Created by Joan Smith
on 2017-4-1.

Given a set of CNV input files for  each cancer type, make a pancan file, including stouffers

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
    opts, args = getopt.getopt(argv[1:], 'hs:i:o:',
    ['help', 'outdir=', 'suffix=', 'indir'])
  except getopt.error, msg:
    usage()

  suffix = None
  indir = None
  outdir = '.'

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-i', '--indir'):
      indir = value
    if option in ('-s', '--suffix'):
      suffix = value

  return outdir, indir, suffix

def make_path(indir, cancer, suffix):
  if os.path.isdir(os.path.join(indir, cancer)):
    return os.path.join(indir, cancer, cancer + suffix)
  else:
    return os.path.join(indir, cancer + suffix)

def make_pancan_df(outdir, indir, suffix):
  pancan = {}
  for cancer in CANCER_TYPES:
    make_path(indir, cancer, suffix)
    path = make_path(indir, cancer, suffix)
    df = pd.read_csv(path, index_col=0, na_values=[' NA']) # sometimes rpy2 gives back this monstrosity of a NaN value.
    if '\'' not in df.index[0]:
      # ADD ' to index
      df = df.reset_index()
      print 'ADD apostrophe'
      df['gene'] = '\'' + df['gene']
      df = df.set_index('gene')
    df_clean = df.drop([u'Chromosome', u'Location'], axis=1)
    pancan[cancer] = df_clean.mean(axis=1).astype(float)

  pancan['Chromosome'] = df['Chromosome']
  pancan['Location'] = df['Location']
  pancan_df = pd.DataFrame(pancan)
  outpath = os.path.join(outdir, 'pancan.csv')
  pancan_df.to_csv(outpath, index_label='gene')

def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, indir, suffix = get_options(argv)
    make_pancan_df(outdir, indir, suffix)

if __name__ == "__main__":
  main()
