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

def get_cna_amplified_ttest(models_treated, other=None):
  print models_treated.name
  duplicate_genes = other.columns[other.columns.duplicated()].values
  if models_treated.name in ['ArmLevelCNScore', 'FocalCNScore']:
    return None
  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')
  ttests = {}
  for c in other:
    if c in duplicate_genes:
      continue
    not_amplified = combined_treatment_result[combined_treatment_result[c] < 0.3]
    amplified = combined_treatment_result[combined_treatment_result[c] >= 0.3]
    if len(not_amplified) >= 7 and len(amplified) >= 7:
      not_amplified_resp  = combined_treatment_result['BestAvgResponse'][not_amplified.index]
      amplified_resp  = combined_treatment_result['BestAvgResponse'][amplified.index]
      if amplified_resp.mean() < 10:
        ttest_res = stats.ttest_ind(not_amplified_resp, amplified_resp, equal_var=False)
        ttests[c] = {'t': ttest_res[0], 'p': ttest_res[1]}

  return pd.DataFrame(ttests)

def get_cna_deleted_ttest(models_treated, other=None):
  duplicate_genes = other.columns[other.columns.duplicated()].values
  print models_treated.name
  if models_treated.name in ['ArmLevelCNScore']:
    return None
  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')
  ttests = {}
  for c in other:
    if c in duplicate_genes:
      continue
    not_deleted = combined_treatment_result[combined_treatment_result[c] > -0.3]
    not_deleted = not_deleted[c]
    deleted = combined_treatment_result[combined_treatment_result[c] <= -0.3]
    if len(not_deleted) >= 7 and len(deleted) >= 7:
      not_deleted_resp  = combined_treatment_result['BestAvgResponse'][not_deleted.index]
      deleted_resp  = combined_treatment_result['BestAvgResponse'][deleted.index]
      if deleted_resp.mean() < 10:
        ttest_res = stats.ttest_ind(not_deleted_resp, deleted_resp, equal_var=False)
        ttests[c] = {'t': ttest_res[0], 'p': ttest_res[1]}

  return pd.DataFrame(ttests)

def get_mutation_ttests(models_treated, other=None):
  print models_treated.name
  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')
  ttests = {}
  for c in other:
    if (combined_treatment_result[c] != 0).sum() >= 7:
      not_mutated = combined_treatment_result[combined_treatment_result[c] == 0].index
      not_mutated_resp = combined_treatment_result['BestAvgResponse'][not_mutated]
      mutated = combined_treatment_result[combined_treatment_result[c] == 1].index
      mutated_resp = combined_treatment_result['BestAvgResponse'][mutated]
      if mutated_resp.mean() < 10:
        ttest_res = stats.ttest_ind(not_mutated_resp, mutated_resp, equal_var=False)
        ttests[c] =  {'t': ttest_res[0], 'p': ttest_res[1]}

  return pd.DataFrame(ttests)

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
  xls = pd.ExcelFile(infile)
  results = xls.parse('PCT curve metrics')

  cna = prep_cna_data(xls)

  cna_corrs = results.groupby(['Treatment']).apply(get_cna_deleted_ttest, other=cna)
  cna_corrs.T.to_csv(os.path.join(outdir, 'xenograft_cna_deleted_lt_10_cutoff_ttests.csv'))

  cna_amplification = results.groupby(['Treatment']).apply(get_cna_amplified_ttest, other=cna)
  cna_amplification.T.to_csv(os.path.join(outdir, 'xenograft_cna_amplification_lt_10_cutoff_ttests.csv'))

  mutations = prep_mutation_data(xls)
  mutation_corrs = results.groupby(['Treatment']).apply(get_mutation_ttests, other=mutations)
  mutation_corrs.T.to_csv(os.path.join(outdir, 'xenograft_mutation_lt_10_cutoff_ttests.csv'))



def main():
  infile, outdir = get_options()
  do_work(infile, outdir)





if __name__ == "__main__":
  main()

