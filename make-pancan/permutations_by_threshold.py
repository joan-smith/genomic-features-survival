#!/usr/bin/env python
# encoding: utf-8
'''
permutations_by_threshhold.py

Created by Joan Smith
on 2017-6-16.

Given a pancan file,
Permute the pancan file n times, then for each row in each permutation, count the number of zscores that pass
the threshhold. This is the permutation histogram.

Then go through the original pancan file and count how many zscores for each gene have a value greater
than (or less than) the provided threshold. Compare this value to the number of times this count appeared after
permutations

Copyright (c) 2017 . All rights reserved.
'''

import pdb
import pandas as pd
import numpy as np
import argparse
import sys
import os
import collections
from matplotlib import pyplot

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib

THRESHOLD = 4
GREATER_THAN = True

def get_options():
  parser = argparse.ArgumentParser(
        description=('Permute pancan data to find out how often a gene has a zscore that passes a threshold'
        'in multiple cancer types'))
  parser.add_argument('-f', action='store', dest='pancan_file', required=True, type=str)
  parser.add_argument('-o', action='store', dest='outdir', required=False, type=str, default='.')
  parser.add_argument('-p', action='store', dest='permutations', required=False, type=int, default=4)
  namespace = parser.parse_args()

  return namespace.pancan_file,  namespace.outdir, namespace.permutations

def shuffle(df, axis=0):
  s = df.copy()
  s = s.apply(np.random.permutation, axis=axis)
  return s


def calculate_permutations(pancan, permutations):
  gt_counts = []
  lt_counts = []
  for i in range(permutations):
    shuffled =  shuffle(pancan)
    permutation_gt_sums = (shuffled > THRESHOLD).sum(axis=1)
    permutation_gt_counts = permutation_gt_sums.value_counts()

    permutation_lt_sums = (shuffled < (THRESHOLD*-1)).sum(axis=1)
    permutation_lt_counts = permutation_lt_sums.value_counts()

    gt_counts.append(permutation_gt_counts)
    lt_counts.append(permutation_lt_counts)

  all_lt_permutations = pd.concat(lt_counts, axis=1)
  total_lt_counts = all_lt_permutations.sum(axis=1, skipna=True)

  all_gt_permutations = pd.concat(gt_counts, axis=1)
  total_gt_counts = all_gt_permutations.sum(axis=1, skipna=True)

  total_counts = pd.concat([total_gt_counts, total_lt_counts], axis=1)
  total_counts.columns = ['Greater than Threshold', 'Less than -Threshold']

  with open('permutation_threshold_counts_%f.csv' % THRESHOLD, 'w') as out:
    out.write('Threshold, %f\n' % THRESHOLD)
    out.write('Permutations, %d\n' % permutations)
    total_counts.to_csv(out, index_label='Cancer Type Count')

def calculate_raw_counts(pancan, permutations):
  greater_than_sums = (pancan > THRESHOLD).sum(axis=1)
  less_than_sums = (pancan < (THRESHOLD*-1)).sum(axis=1)

  sums = pd.concat([greater_than_sums, less_than_sums], axis=1)
  sums.columns = ['Greater than Threshold', 'Less than -Threshold']

  with open('pancan_threshold_counts_%f.csv' % THRESHOLD, 'w') as out:
    out.write('Threshold, %f\n' % THRESHOLD)
    out.write('Permutations, %d\n\n' % permutations)
    sums.to_csv(out, index_label='Gene')

def main(argv=None):
  pancan_file, outdir, permutations = get_options()
  pancan = pd.read_csv(pancan_file, header=0, index_col=0)
  pancan = pancan.astype(float)
  pancan = pancan.drop('stouffer unweighted', axis=1)

  calculate_raw_counts(pancan, permutations)
  calculate_permutations(pancan, permutations)


if __name__ == "__main__":
  main()
