#!/usr/bin/env python
# encoding: utf-8
'''

zscores_for_copy_number_when_mutated.py


Created by Joan Smith
on 2017-7-7.

Given a set of interesting genes, get the correlations between copy number and rnaseq.

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

sys.path.append('../../common/')
import utilities as util
import analysis
import mutation_base

def get_options():
  parser = argparse.ArgumentParser(description='Get cnv, rna, and clinical directories. optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_directory')
  parser.add_argument('-r', action='store', dest='rna')
  parser.add_argument('-m', action='store', dest='mutation_dir')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-f', action='store', dest='interesting_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.cnv_directory, namespace.rna, namespace.mutation_dir,
          namespace.clinical_directory, namespace.output_directory,
          namespace.interesting_file)


def make_corrs(copy_number, rnaseq, mutation, clinical, outdir, genes):
  cancer_type = util.get_cancer_type(copy_number)

  clinical_data = util.get_clinical_data(clinical)
  cnv = pd.read_csv(copy_number, index_col=0)
  cnv_by_patient = cnv.transpose()


  rnaseq =  pd.read_csv(rnaseq, low_memory=False, sep='\t')
  rnaseq = rnaseq.drop([0])
  rnaseq = rnaseq.set_index('Hybridization REF').astype(np.float)
  rnaseq = rnaseq.transpose().reset_index()
  rnaseq = util.maybe_clear_non_01s(rnaseq, 'index', cancer_type)
  rnaseq = util.add_identifier_column(rnaseq, 'index')
  rnaseq_clean = rnaseq.set_index('identifier').drop('index', 1).astype(np.float)
  rnaseq_log2 = rnaseq_clean.apply(np.log2)
  rnaseq_clipped_log2 = np.clip(rnaseq_log2, 0, np.inf)
  rna_cnv = cnv_by_patient[genes['Gene']].join(rnaseq_clipped_log2, how='inner')

  mutation = mutation_base.prep_mutation_data(mutation, clinical_data)
  print mutation.index

  included_patients = set(list(mutation.index)) & set(list(rna_cnv.index))

  rna_cnv = rna_cnv.loc[included_patients]

  rna_cnv.T.to_csv(os.path.join(outdir, cancer_type + '_cnv_rnaseq_data.csv'))

  corr_dict = {}
  for gene in genes['Gene']:
    corr = rna_cnv.corrwith(rna_cnv[gene]).drop(genes['Gene'])
    corr_dict[cancer_type + '_' + gene] = corr

  return pd.DataFrame(corr_dict)

def multiprocess_data(args):
  copy_number = args[0]
  rnaseq = args[1]
  mutation = args[2]
  clinical = args[3]
  outdir = args[4]
  genes = args[5]
  corr_df = make_corrs(copy_number, rnaseq, mutation, clinical, outdir, genes)
  print corr_df
  return corr_df


def main(argv=None):
  cnv_dir, rna, mutation_dir, clinical_dir, outdir, input_file = get_options()
  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  interesting_genes = pd.read_csv(input_file, comment='#')
  interesting_genes['Gene'] = '\'' + interesting_genes['Gene']

  zscore_inputs = []
  corr_results = []
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    cancer_type_genes = interesting_genes[interesting_genes['Cancer Type'] == cancer_type]
    if len(cancer_type_genes) == 0:
      continue
    print cancer_type
    print cancer_type_genes
    rnaseq = glob.glob(os.path.join(rna, cancer_type + '*'))[0]
    mutation = glob.glob(os.path.join(mutation_dir, cancer_type + '*'))[0]
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]


    zscore_inputs.append([cnv, rnaseq, mutation, clinical, outdir, cancer_type_genes])
    # corr_results.append(multiprocess_data([cnv, rnaseq, mutation, clinical, outdir, cancer_type_genes]))

  p = Pool(4)
  corr_results = p.map(multiprocess_data, zscore_inputs)
  df = pd.concat(corr_results, verify_integrity=True, axis=1)
  print df
  df.to_csv(os.path.join(outdir, 'corr_results.csv'))






if __name__ == "__main__":
  main()
