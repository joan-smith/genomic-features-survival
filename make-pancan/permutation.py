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
  s = s.apply(np.random.permutation, axis=axis)
  return s

def calculate(df):
  np.random.seed()
  shuffled = shuffle(df)
  stouffer = analysis.stouffer_unweighted(shuffled)
  return stouffer.values

def do_shuffle_work(pancan_file, outdir, permutations, threads, show):
  df = pd.read_csv(pancan_file, index_col=0)
  histogram_data = []
  pool = multiprocessing.Pool(threads)
  histogram_data = pool.map(calculate, [df]*permutations)
  histogram_data = np.array(histogram_data).flatten()
  util.make_histogram(histogram_data, outdir, show)

def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, infile, permutations, threads, show = get_options(argv)
    do_shuffle_work(infile, outdir, permutations, threads, show)


if  __name__ == '__main__':
  main()
