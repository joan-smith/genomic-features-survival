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


def get_options():
  parser = argparse.ArgumentParser(description='Get mutation, cnv, and clinical directories. optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_directory')
  parser.add_argument('-c', action='store', dest='raw_clinical_directory')
  parser.add_argument('-f', action='store', dest='interesting_genes_file')
  parser.add_argument('-s', action='store', dest='structural_breaks')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.cnv_directory, ns.structural_breaks, ns.interesting_genes_file, ns.output_directory)

def make_cn_zscores(copy_number, breaks, interesting_genes=None, outdir='.'):
  clinical_data = pd.read_csv(breaks, index_col=0)
  cnv = pd.read_csv(copy_number, index_col=0)
  cnv_by_patient = cnv.transpose()
  clinical_and_cnv = cnv_by_patient.join(clinical_data, how='inner')

  cancer_type = util.get_cancer_type(copy_number)
  print cancer_type

  if interesting_genes is not None:
    relevant_genes = ('\'' + interesting_genes[interesting_genes['Cancer Type'] == cancer_type]['Gene']).values
    relevant_genes = list(relevant_genes) + ['breaks', 'time', 'censor']
    print relevant_genes
    clinical_and_cnv = clinical_and_cnv[relevant_genes]
    clinical_and_cnv.to_csv(os.path.join(outdir, cancer_type + '_interesting_genes_data.csv'))
    return


  results = []
  for gene in clinical_and_cnv:
    if gene in ('time', 'censor'): # skip metadata
      continue
    if clinical_and_cnv[gene].count() > 10:
      cox_dict = analysis.do_multivariate_cox(clinical_and_cnv.time,
                                              clinical_and_cnv.censor,
                                              clinical_and_cnv[gene],
                                              clinical_and_cnv[['breaks']])
      cox_dict['gene'] = gene
      results.append(cox_dict)
  results_df = pd.DataFrame(results)
  results_df = results_df.set_index('gene')
  results_df.to_csv(os.path.join(outdir, cancer_type + '_cn_cox.csv'))

def multiprocess_cn_zscores(args):
  copy_number, breaks, interesting_genes, outdir = args
  make_cn_zscores(copy_number, breaks, outdir)

def main(argv=None):
  cnv_dir, structural_breaks, interesting_genes_file, outdir = get_options()
  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  interesting_genes = None
  if interesting_genes_file:
    interesting_genes = pd.read_csv(interesting_genes_file)


  zscore_inputs = []
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    breaks = glob.glob(os.path.join(structural_breaks, '*' + cancer_type + '*'))[0]
    zscore_inputs.append([cnv, breaks, interesting_genes, outdir])
    make_cn_zscores(cnv, breaks, interesting_genes, outdir)

  p = Pool(4)
  p.map(multiprocess_cn_zscores, zscore_inputs)

if __name__ == "__main__":
  main()
