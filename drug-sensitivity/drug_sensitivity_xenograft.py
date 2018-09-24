#!/usr/bin/env python
# encoding: utf-8
'''
drug_sensitivity.py

Created by Joan Smith
on 2017-10-21

Given experimental data, calculate pearson correlations between gene expression and best avg response

Copyright (c) 2018. All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np
import scipy.stats as stats
from IPython import embed

INCLUDED_MUTATIONS = ['MutNovel', 'MutKnownFunctional', 'MutLikelyFunctional']

def get_options():
  parser = argparse.ArgumentParser(description='Produce intermediate file for analysis.')
  parser.add_argument('-i', action='store', dest='input_file', required=True)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_file, namespace.output_directory)

def get_correlation(models_treated, other=None):
  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')
  return combined_treatment_result.corrwith(combined_treatment_result['BestAvgResponse'])

def prep_cna_data(xls):
  cna = xls.parse('copy number', index_col=0).T
  cna = cna.divide(2)
  cna = cna.apply(np.log2)
  return cna

def prep_mutation_data(xls):
  mutation = xls.parse('pdxe_mut_and_cn2')
  mutation = mutation[mutation['Category'].isin(INCLUDED_MUTATIONS)]
  gene_mutation_df = mutation.groupby(['Gene']).apply(mutations_for_gene)
  gene_mutation_df.index.set_names(['Gene', 'Model'], inplace=True)
  gene_mutation_df = gene_mutation_df.reset_index()
  gene_patient_mutations = gene_mutation_df.pivot(index='Gene', columns='Model', values='mutated')
  return gene_patient_mutations.fillna(0).T

def mutations_for_gene(df):
  mutated_patients = df['Sample'].unique()
  return pd.DataFrame({'mutated': np.ones(len(mutated_patients))}, index=mutated_patients)


def do_work(infile, outdir):
  print infile
  xls = pd.ExcelFile(infile)
  results = xls.parse('PCT curve metrics')

  # cna = prep_cna_data(xls)
  # cna_corrs = results.groupby(['Treatment']).apply(get_correlation, other=cna)
  # cna_corrs.T.to_csv(os.path.join(outdir, 'CNA_correlations.csv'))
  #
  # gene_expr = xls.parse('RNAseq_fpkm', index_col=0).T
  # gene_expr_corrs = results.groupby(['Treatment']).apply(get_correlation, other=gene_expr)
  # gene_expr_corrs.T.to_csv(os.path.join(outdir, 'gene_expr_correlations.csv'))

  mutations = prep_mutation_data(xls)
  mutations.to_csv('raw_mutations_data.csv')
  embed()
  mutation_corrs = results.groupby(['Treatment']).apply(get_correlation, other=mutations)
  mutation_corrs.T.to_csv(os.path.join(outdir, 'mutation_correlations.csv'))



def main():
  infile, outdir = get_options()
  do_work(infile, outdir)





if __name__ == "__main__":
  main()

