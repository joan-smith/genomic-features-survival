#!/usr/bin/env python
# encoding: utf-8
'''
make pancan.py

Created by Joan Smith
on 2017-4-1.

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

CANCER_TYPES_ABRIDGED =  ['BLCA', 'BRCA',  'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'LUAD',
                'LUSC', 'OV', 'PAAD', 'PRAD', 'SARC', 'STES']

ICGC_TYPES = ['COCA-CN', 'EOPC-DE', 'ESCA-CN', 'GACA-CN', 'LICA-FR', 'LINC-JP', 'LIRI-JP',
              'MELA-AU', 'ORCA-IN', 'OV-AU', 'PACA-AU',
             'PACA-AU', 'PACA-CA', 'PBCA-DE', 'PRAD-UK', 'RECA-EU', 'SKCA-BR']

CBIOPORTAL_TYPES = ['BLCA', 'BRCA', 'LIHC', 'LUAD', 'PRAD']


def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'hs:i:o:mgc',
    ['help', 'outdir=', 'suffix=', 'indir=', 'metagene', 'icgc', 'cbioportal'])
  except getopt.error, msg:
    usage()

  suffix = None
  indir = None
  metagene = False
  outdir = None
  icgc = False
  cbioportal = False

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-i', '--indir'):
      indir = value
    if option in ('-s', '--suffix'):
      suffix = value
    if option in ('-m', '--metagene'):
      metagene = True
    if option in ('-g', '--icgc'):
      icgc = True
    if option in ('-c', '--cbioportal'):
      cbioportal = True

  if not outdir:
    outdir = indir

  return outdir, indir, suffix, metagene, icgc, cbioportal

def make_path(indir, cancer, suffix):
  if os.path.isdir(os.path.join(indir, cancer)):
    return os.path.join(indir, cancer, cancer + suffix)
  else:
    return os.path.join(indir, cancer + suffix)

def make_pancan_df(outdir, indir, suffix, metagene=False, icgc=False, cbioportal=False):
  pancan = {}
  zscore_header = 'zscore'
  outfile = 'pancan.csv'
  if metagene:
    print 'Calculating for metagene'
    zscore_header = 'metagene-zscore'
    outfile = 'metagene_pancan.csv'

  cancer_types = CANCER_TYPES
  if icgc:
    cancer_types = ICGC_TYPES
  elif cbioportal:
    cancer_types = CBIOPORTAL_TYPES

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
      pancan[cancer] = df[zscore_header].astype(float)

  pancan_df = pd.DataFrame(pancan)
  pancan_df['stouffer unweighted'] = analysis.stouffer_unweighted(pancan_df)
  outpath = os.path.join(outdir, outfile)
  pancan_df.to_csv(outpath, index_label='gene')

def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, indir, suffix, metagene, icgc, cbioportal = get_options(argv)
    make_pancan_df(outdir, indir, suffix, metagene, icgc, cbioportal)

if __name__ == "__main__":
  main()
