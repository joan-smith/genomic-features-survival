#!/usr/bin/env python
# encoding: utf-8
'''
permutation.py

Created by Joan Smith
on 2016-02-20.

Given a pancan file, permute the results and make a histogram

Copyright (c) 2017 . All rights reserved.
'''

import getopt
import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot
import multiprocessing

sys.path.append('../common/')
import analysis

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'hi:o:p:t:s',
    ['help', 'outdir=', 'suffix=', 'infile=', 'permutations=', 'threads=', 'show'])
  except getopt.error, msg:
    usage()

  infile = None
  outdir = '.'
  permutations = 4
  threads = 4
  show = False

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-i', '--infile'):
      infile = value
    if option in ('-p', '--permutations'):
      permutations = int(value)
    if option in ('-t', '--threads'):
      threads = int(value)
    if option in ('-s', '--show'):
      show = True

  return outdir, infile, permutations, threads, show

def shuffle(df, axis=0):
  s = df.copy()
  s = s.apply(np.random.permutation, axis=0)
  return s

def calculate(df):
  shuffled = shuffle(df, axis=1)
  stouffer = analysis.stouffer_unweighted(shuffled)
  return stouffer.values

def do_shuffle_work(pancan_file, outdir, permutations, threads, show):
  df = pd.read_csv(pancan_file, index_col=0)
  df = df.drop(u'stouffer unweighted', axis=1) # clear the initial stouffer column

  histogram_data = []
  pool = multiprocessing.Pool(threads)
  histogram_data = pool.map(calculate, [df]*permutations)
  histogram_data = np.array(histogram_data).flatten()
  minimum = min(histogram_data)
  maximum = max(histogram_data)

  pyplot.title(outdir)
  hist_counts, bins, _ = pyplot.hist(histogram_data, bins=100)
  hist_counts = np.append(hist_counts, np.nan)
  bucketed_df = pd.DataFrame({'bins': bins, 'count': hist_counts})

  with open(os.path.join(outdir, 'panplatform_histogram_data.csv'), 'w') as out:
    out.write('permutations,' + str(permutations) + '\n')
    out.write('min,' + str(minimum) + '\n')
    out.write('max,' + str(maximum) + '\n')
    bucketed_df.to_csv(out, index=False)

  pyplot.savefig(os.path.join(outdir, 'panplatform_histogram.png'))

  print outdir
  print minimum
  print maximum
  if show:
    pyplot.show()

def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, infile, permutations, threads, show = get_options(argv)
    do_shuffle_work(infile, outdir, permutations, threads, show)


if  __name__ == '__main__':
  main()
