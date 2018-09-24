#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-08.

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

sys.path.append('../mutation-analysis')
import zscores_for_mutants as mutation_zscores

MUTATION_PERCENT = 0.02

def get_options():
  parser = argparse.ArgumentParser(description='Get mutation, and clinical directories. optional output dir')
  parser.add_argument('-m', action='store', dest='mutation_directory')
  parser.add_argument('-c', action='store', dest='raw_clinical_directory')
  parser.add_argument('-s', action='store', dest='structural_breaks')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.mutation_directory, ns.raw_clinical_directory, ns.structural_breaks,
    ns.output_directory)

def make_zscores(mutation, clinical, breaks, outdir):
  clinical_data = util.get_clinical_data(clinical)
  mut = mutation_base.prep_mutation_data(mutation, clinical_data)

  cancer_type = util.get_cancer_type(mutation)
  print cancer_type

  structural_breaks = pd.read_csv(breaks, index_col=0)
  structural_breaks = structural_breaks.astype(int)
  mut_and_breaks = mut.join(structural_breaks, how='inner')
  num_patients = len(mut_and_breaks)

  results = []
  for gene in mut_and_breaks:
    if gene in ('time', 'censor', 'breaks'): # skip metadata
      continue
    num_mutations = mut_and_breaks[gene].sum()
    if num_mutations >= MUTATION_PERCENT * num_patients:
      cox_dict = analysis.do_multivariate_cox(mut_and_breaks.time,
                                              mut_and_breaks.censor,
                                              mut_and_breaks[gene],
                                              mut_and_breaks[['breaks']])
      cox_dict['gene'] = gene
      results.append(cox_dict)
  results_df = pd.DataFrame(results)
  results_df = results_df.set_index('gene')
  results_df.to_csv(os.path.join(outdir, cancer_type + '_mut_cox.csv'))

def multiprocess_zscores(args):
  mut, clinical, breaks, outdir = args
  make_zscores(mut, clinical, breaks, outdir)

def main(argv=None):
  mutation_dir, clinical_dir, structural_breaks, outdir = get_options()
  mut_files = os.listdir(mutation_dir)
  mut_files = util.remove_extraneous_files(mut_files)
  mut_files = [os.path.join(mutation_dir, i) for i in mut_files]

  zscore_inputs = []
  for mut in mut_files:
    if '_' in mut:
      continue
    cancer_type = util.get_cancer_type(mut)
    print cancer_type
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]
    breaks = glob.glob(os.path.join(structural_breaks, '*' + cancer_type + '*'))[0]
    zscore_inputs.append([mut, clinical, breaks, outdir])
    # make_zscores(mut, clinical, breaks, outdir)

  p = Pool(4)
  p.map(multiprocess_zscores, zscore_inputs)

if __name__ == "__main__":
  main()
