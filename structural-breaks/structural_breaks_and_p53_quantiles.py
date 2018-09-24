#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-07.

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
from scipy import stats


sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base

def get_options():
  parser = argparse.ArgumentParser(description='Make quantiles of structural breaks data for p53 mutation statuses')
  parser.add_argument('-m', action='store', dest='mutation_dir')
  parser.add_argument('-c', action='store', dest='structural_breaks_clinical')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.mutation_dir, ns.structural_breaks_clinical, ns.output_directory)

def main():
  mutation_dir, clinical_dir, outdir = get_options()
  mutation_files = os.listdir(mutation_dir)
  mutation_files = util.remove_extraneous_files(mutation_files)

  results = pd.DataFrame()
  for mut in mutation_files:
    if '_' in mut:
      continue
    cancer_type = util.get_cancer_type(mut)
    print cancer_type
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]

    clinical_data = pd.read_csv(clinical, index_col=0)
    mutation = mutation_base.prep_mutation_data(os.path.join(mutation_dir, mut), clinical_data)
    data = mutation[['\'TP53']].join(clinical_data, how='inner')
    print data

    wt_as = data[data['\'TP53'] == 0]['breaks']
    mut_as = data[data['\'TP53'] != 0]['breaks']

    wt_q = wt_as.quantile([0.10, 0.25, 0.50, 0.75, 0.90])
    mut_q = mut_as.quantile([0.10, 0.25, 0.50, 0.75, 0.90])

    statistic, p = stats.mannwhitneyu(wt_as, mut_as)

    wt_q['cancer_type'] = cancer_type
    wt_q['mut?'] = 'wt'
    mut_q['cancer_type'] = cancer_type
    mut_q['mut?'] = 'mut'
    wt_q['mann-whitney-p'] = p

    results = results.append(wt_q)
    results = results.append(mut_q)

  results = results.set_index(['cancer_type', 'mut?'])
  results.to_csv(os.path.join(outdir, 'breaks_and_p53_quantiles.csv'))


if __name__ == "__main__":
  main()
