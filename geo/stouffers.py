#!/usr/bin/env python
# encoding: utf-8
'''
stouffers.py

Created by Joan Smith
on 2017-7-6.

Given a geo pancan file, calculate stouffers for each type, and then stouffers on that
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os

sys.path.append('../common/')
import analysis

ACROSS_CANCER_TYPE_CUTOFF = 0.5
WITHIN_CANCER_TYPE_CUTOFF = 0.5


def get_options(argv):
  parser = argparse.ArgumentParser(description='Produce a stouffers file')
  parser.add_argument('-i', action='store', dest='input_file')
  parser.add_argument('-o', action='store', dest='outdir', default=None)
  namespace = parser.parse_args()

  return (namespace.input_file, namespace.outdir)


def cancertype(ind):
  return ind[2]

def make_pancan_df(outdir, indir, suffix, metagene=False, icgc=False, cbioportal=False):
  pancan_df = pd.DataFrame(pancan)
  pancan_df['stouffer unweighted'] = analysis.stouffer_unweighted(pancan_df)
  outpath = os.path.join(outdir, outfile)
  pancan_df.to_csv(outpath, index_label='gene')

def make_stouffer(input_file, outdir):
  df = pd.read_csv(input_file, header=None, index_col=0, low_memory=False)
  df.columns = pd.MultiIndex.from_arrays([df.loc['GPL'], df.loc['GSE'], df.loc['type'], df.loc['patient count']])
  df = df.drop(['GPL', 'GSE', 'type', 'patient count'])
  total_gses = df.shape[1]
  df_cutoff = df[df.count(axis=1) >= total_gses*ACROSS_CANCER_TYPE_CUTOFF]

  cancer_type_groups = df_cutoff.groupby(by=cancertype, axis=1)
  ctype_stouffer = {}
  for t, g in cancer_type_groups:
    print t, g.shape
    gses_in_type = g.shape[1]
    g = g[g.count(axis=1) >= gses_in_type*WITHIN_CANCER_TYPE_CUTOFF]
    print g.shape

    ctype_stouffer[t] = analysis.stouffer_unweighted(g.astype(float))

  ctype_stouffer_df = pd.DataFrame(ctype_stouffer)
  ctype_stouffer_df['cross-type stouffer'] = analysis.stouffer_unweighted(ctype_stouffer_df)
  ctype_stouffer_df.to_csv(os.path.join(outdir, 'cancer_type_stouffers.csv'), index_label='Cancer Type')


def main(argv=None):
  input_file, outdir = get_options(argv)
  if outdir == None:
    outdir = os.path.dirname(input_file)

  make_stouffer(input_file, outdir)


if __name__ == "__main__":
  main()
