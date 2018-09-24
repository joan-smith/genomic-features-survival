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


def get_options():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('-i', action='store', dest='input_directory', required=True)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.output_directory)

def get_cna_amplified_ttest(models_treated, other=None):
  print models_treated.name
  combined_treatment_result = models_treated.T.join(other, how='inner')
  ttests = {}
  for c in other:
    not_amplified = combined_treatment_result[combined_treatment_result[c] < 3]
    amplified = combined_treatment_result[combined_treatment_result[c] >= 3]
    if len(not_amplified) > 4 and len(amplified) > 4:
      not_amplified_resp  = combined_treatment_result[models_treated.name][not_amplified.index]
      amplified_resp  = combined_treatment_result[models_treated.name][amplified.index]
      ttest_res = stats.ttest_ind(not_amplified_resp, amplified_resp, equal_var=False)
      ttests[c] = {'t': ttest_res[0], 'p': ttest_res[1]}

  return pd.DataFrame(ttests)

def get_cna_deleted_ttest(models_treated, other=None):
  print models_treated.name
  combined_treatment_result = models_treated.T.join(other, how='inner')
  ttests = {}
  for c in other:
    not_deleted = combined_treatment_result[combined_treatment_result[c].isin([0,1])]
    deleted = combined_treatment_result[~combined_treatment_result[c].isin([0,1])]
    if len(not_deleted) > 4 and len(deleted) > 4:
      not_deleted_resp  = combined_treatment_result[models_treated.name][not_deleted.index]
      deleted_resp  = combined_treatment_result[models_treated.name][deleted.index]
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
    if (combined_treatment_result[c] != 0).sum() >= 4:
      not_mutated = combined_treatment_result[combined_treatment_result[c] == 0].index
      not_mutated_resp = combined_treatment_result['BestAvgResponse'][not_mutated]
      mutated = combined_treatment_result[combined_treatment_result[c] == 1].index
      mutated_resp = combined_treatment_result['BestAvgResponse'][mutated]
      if mutated_resp.mean() < 20:
        ttest_res = stats.ttest_ind(not_mutated_resp, mutated_resp, equal_var=False)
        ttests[c] =  {'t': ttest_res[0], 'p': ttest_res[1]}

  return pd.DataFrame(ttests)



def get_copy_number(gene):
  copy_number = gene.str.split(',', expand=True)[0]
  copy_number = copy_number.astype(int)
  copy_number = copy_number.replace(-1, np.nan)
  return copy_number

def prep_cnas(indir):
  cnas_xls = pd.ExcelFile(os.path.join(indir, 'Gene_level_CN.xlsx'))
  cnas = cnas_xls.parse('Gene_level_CN')
  cnas.columns = list(cnas.columns[:4]) + list(cnas.loc[0][4:])
  cnas = cnas.drop(0)
  cnas = cnas.drop(['start', 'stop', 'chr'], axis=1)
  cnas.columns = ['gene'] + list(cnas.columns[1:].astype(int).astype(str))

  cnas = cnas.set_index('gene').T

  cnas = cnas.apply(get_copy_number)
  return cnas

def prep_annotation(indir):
  annotation_xls = pd.ExcelFile(os.path.join(indir, 'TableS1E - annotations.xlsx'))
  annotation = annotation_xls.parse('TableS1E-CellLines', parse_cols='C,K', header=1)
  annotation.columns = ['COSMIC', 'cancer_type']
  annotation = annotation.dropna(how='all')
  # COSMIC is being parsed as a float, but it needs to be a string.
  annotation['COSMIC'] = annotation['COSMIC'].astype(int).astype(str)
  annotation = annotation.set_index('COSMIC')
  return annotation

def run_cna(treatments, indir, outdir):
  annotation = prep_annotation(indir)
  cnas = prep_cnas(indir)

  # do pancan analsyes
  pancan_amplified = treatments.groupby(level=0).apply(get_cna_amplified_ttest,
                                  other=cnas)
  pancan_amplified.T.to_csv(os.path.join(outdir, 'cell_line_amplification_cna_ttest.csv'))

  pancan_deleted = treatments.groupby(level=0).apply(get_deletedd_ttest,
                                  other=cnas)
  pancan_deleted.T.to_csv(os.path.join(outdir, 'cell_line_deleted__cna_ttest.csv'))


  # by cancer type analyses
  cancer_type_groups = cnas.join(annotation).groupby('cancer_type')
  for cancer_type, cnas_of_type in cancer_type_groups:
    cancer_type = cancer_type.replace('/', '')
    cnas_of_type = cnas_of_type.drop('cancer_type', axis=1)
    print cancer_type, cnas_of_type.shape[0]
    if cnas_of_type.shape[0] >= 10:
      amplified_cnas = treatments.groupby(level=0).apply(get_cna_amplified_ttest,
                                               other=cnas_of_type)
      amplified_cnas.T.to_csv(os.path.join(outdir, cancer_type + '_cell_line_amplified_cna_ttest.csv'))

      deleted_cnas = treatments.groupby(level=0).apply(get_cna_deleted_ttest,
                                               other=cnas_of_type)
      deleted_cnas.T.to_csv(os.path.join(outdir, cancer_type + '_cell_line_deleted_cna_ttest.csv'))


def clean_treatment(t):
  t = t.drop('COSMIC', axis=1)
  return t.mean()


def prep_treatment(indir):
  treatment_xls = pd.ExcelFile(os.path.join(indir, 'drug sensitivity copy for ingestion.xlsx'))
  treatments = treatment_xls.parse('TableS4A-IC50s', index_col=0, header=None)
  treatments = treatments.T
  t = treatments.groupby(['COSMIC']).apply(clean_treatment)
  t.columns = t.columns.astype(int).astype(str)
  return t


def main():
  indir, outdir = get_options()

  treatments = prep_treatment(indir)
  run_cna(treatments, indir, outdir)


if __name__ == "__main__":
  main()

