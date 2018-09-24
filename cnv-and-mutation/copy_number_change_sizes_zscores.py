#!/usr/bin/env python
# encoding: utf-8
'''

Created by Joan Smith
on 2018-2-12.

Given the copy-number-change sizes files generated from the sister script,
calculate z-scores first excluding broad mutations, then excluding focal mutations.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import collections
import glob
import rpy2
from multiprocessing import Pool

sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base

FOCAL_CUTOFF = 3e6 # any change larger than 3million bps is a broad change.

def get_options():
  parser = argparse.ArgumentParser(description='CN Change size zscores')
  parser.add_argument('-i', action='store', dest='cnv_change_size_dir')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.cnv_change_size_dir, namespace.clinical_directory,
          namespace.output_directory)


def calculate_zscores(input_file):
  input_data = pd.read_csv(input_file, index_col=0)
  input_data = input_data.dropna(subset=['time', 'censor', 'copy number'], how='any')

  input_data['focal'] = input_data.continuous_len <= FOCAL_CUTOFF
  input_data['broad'] = input_data.continuous_len > FOCAL_CUTOFF
  # print input_data
  focal_zscore = analysis.do_cox(input_data.time, input_data.censor,
                                            input_data['focal'])
  focal_zscore['focal_count'] = input_data.focal.sum()
  print focal_zscore

  return focal_zscore


def calculate_broad_change_zscores(input_file):
  input_data = pd.read_csv(input_file, index_col=0)
  input_data = input_data.dropna(subset=['time', 'censor', 'copy number'], how='any')

  input_data['broad'] = input_data.continuous_len > FOCAL_CUTOFF
  # print input_data
  broad_zscore = analysis.do_cox(input_data.time, input_data.censor,
                                            input_data['broad'])
  broad_zscore['broad_count'] = input_data.broad.sum()
  print broad_zscore

  return broad_zscore

def calculate_broad_change_restricted_zscores(input_file):
  input_data = pd.read_csv(input_file, index_col=0)
  input_data = input_data.dropna(subset=['time', 'censor', 'copy number'], how='any')
  print input_data.shape

  # ignore patients that have a focal change
  input_data = input_data.drop(input_data[input_data.continuous_len <= FOCAL_CUTOFF].index)

  input_data['broad'] = input_data.continuous_len > FOCAL_CUTOFF
  broad_restricted_zscore = analysis.do_cox(input_data.time,
                                input_data.censor, input_data['broad'])
  broad_restricted_zscore['broad_count'] = input_data.broad.sum()
  print broad_restricted_zscore
  return broad_restricted_zscore


def calculate_focal_change_restricted_zscores(input_file):
  input_data = pd.read_csv(input_file, index_col=0)
  input_data = input_data.dropna(subset=['time', 'censor', 'copy number'], how='any')
  print input_data.shape

  # ignore patients that have a broad change
  input_data = input_data.drop(input_data[input_data.continuous_len > FOCAL_CUTOFF].index)
  print input_data.shape

  input_data['focal'] = input_data.continuous_len <= FOCAL_CUTOFF
  focal_restricted_zscore = analysis.do_cox(input_data.time,
                                input_data.censor, input_data['focal'])
  focal_restricted_zscore['focal_count'] = input_data.focal.sum()
  print focal_restricted_zscore
  return focal_restricted_zscore

def calculate_any_change_zscores(input_file):
  input_data = pd.read_csv(input_file, index_col=0)
  input_data = input_data.dropna(subset=['time', 'censor', 'copy number'], how='any')
  print input_data.shape

  input_data['any_change'] = ~np.isnan(input_data.continuous_len)
  any_change_zscore = analysis.do_cox(input_data.time,
                                input_data.censor, input_data['any_change'])
  any_change_zscore['any_change_count'] = input_data.any_change.sum()
  print any_change_zscore
  return any_change_zscore


def unused():
  # exclude_broad = input_data[(input_data.continuous_len <= FOCAL_CUTOFF) |
  #                            np.isnan(input_data.continuous_len)]
  # exclude_broad_zscore = analysis.do_cox(exclude_broad.time,
  #                                       exclude_broad.censor, exclude_broad['copy number'])
  # exclude_broad_zscore['included focal count'] = (
  #           exclude_broad.continuous_len < FOCAL_CUTOFF).sum()
  #
  #
  # exclude_focal = input_data[(input_data.continuous_len > FOCAL_CUTOFF) |
  #                            np.isnan(input_data.continuous_len)]
  # exclude_focal_zscore = analysis.do_cox(exclude_focal.time,
  #                                       exclude_focal.censor, exclude_focal['copy number'])
  # exclude_focal_zscore['included broad count'] = (
  #           exclude_focal.continuous_len >= FOCAL_CUTOFF).sum()
  return {'exclude_broad': exclude_broad_zscore, 'exclude_focal': exclude_focal_zscore}



def multiprocess_zscores(args):
  input_file = args[0]
  cancer_type = args[1]
  gene = args[2]
  zscores = calculate_any_change_zscores(input_file)
  return {cancer_type + '_' + gene: zscores}

def main(argv=None):
  cn_change_size_dir, clinical_dir, outdir = get_options()
  input_files = os.listdir(cn_change_size_dir)
  input_files = util.remove_extraneous_files(input_files)
  input_files = [os.path.join(cn_change_size_dir, i) for i in input_files]

  zscore_inputs = []
  results = []
  for input_file in input_files:
    cancer_type = os.path.split(input_file)[1].split('_')[0]
    gene = os.path.split(input_file)[1].split('_')[1].split('.')[0]
    print cancer_type, gene


    # zscore_inputs.append([input_file, cancer_type, gene])
    results.append(multiprocess_zscores([input_file, cancer_type, gene]))

  #p = Pool(4)
  #results = p.map(multiprocess_zscores, zscore_inputs)
  with open(os.path.join(outdir, 'cox_any_change_results.csv'), 'w') as out:
    formatstr = '{},{},{},{}\n'
    out.write('Cancer Type,Gene,Z Score,Count\n')
    for cox_dict in results:
      cancer_type_gene = cox_dict.keys()[0]
      print cancer_type_gene
      print cox_dict[cancer_type_gene]
      d = cox_dict[cancer_type_gene]
      out.write(formatstr.format(cancer_type_gene.split('_')[0], cancer_type_gene.split('_')[1],
              d['z'], d['any_change_count']))

if __name__ == "__main__":
  main()
