#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-10.

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

def get_options():
  parser = argparse.ArgumentParser(description='Get mutation, cnv, and clinical directories. optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_directory')
  parser.add_argument('-c', action='store', dest='clinical')
  parser.add_argument('-f', action='store', dest='interesting_genes_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.cnv_directory, ns.clinical, ns.interesting_genes_file, ns.output_directory)

def make_cn_zscores(copy_number, clinical, interesting_genes=None, outdir='.'):
  clinical_data = util.get_clinical_data(clinical)
  cnv = pd.read_csv(copy_number, index_col=0)
  cnv_by_patient = cnv.transpose()

  cancer_type = util.get_cancer_type(copy_number)

  relevant_genes = '\'' + interesting_genes.index
  relevant_genes = list(relevant_genes)
  cnv = cnv_by_patient[relevant_genes]

  cnv = cnv.join(clinical_data, how='inner')

  results = []
  for gene in cnv:
    if gene in ('time', 'censor'): # skip metadata
      continue
    if cnv[gene].count() > 10:
      cnv[gene + '_split'] = np.nan
      cnv.loc[cnv[gene] <= -0.3, gene + '_split'] = -1
      cnv.loc[cnv[gene].between(-0.3, 0.3), gene + '_split'] = 0
      cnv.loc[cnv[gene] >= 0.3, gene + '_split'] = 1

      cox_dict = analysis.do_cox(cnv.time,
                                cnv.censor,
                                cnv[gene + '_split'])
      cox_dict['gene'] = gene
      cox_dict['cancer_type'] = cancer_type
      results.append(cox_dict)
  cnv.to_csv(os.path.join(outdir, cancer_type + '_trichotomized.csv'))
  return results

def main(argv=None):
  cnv_dir, clinical, interesting_genes_file, outdir = get_options()
  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  interesting_genes = pd.read_csv(interesting_genes_file, index_col=0, header=None)

  results = []
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    clinical_file = glob.glob(os.path.join(clinical, '*' + cancer_type + '*'))[0]
    results += make_cn_zscores(cnv, clinical_file, interesting_genes, outdir)

  results_df = pd.DataFrame(results)
  results_df = results_df.set_index(['cancer_type', 'gene'])
  results_df.to_csv(os.path.join(outdir, 'trichotomized_copy_number_zscores.csv'))



if __name__ == "__main__":
  main()
