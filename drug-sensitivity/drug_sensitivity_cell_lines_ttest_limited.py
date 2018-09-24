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
  parser.add_argument('-l', action='store', dest='limited_query', required=True)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.limited_query, namespace.output_directory)

def get_cna_amplified_ttest(models_treated, other=None):
  combined_treatment_result = models_treated.T.join(other, how='inner')

  if len(other.columns) != 1:
    print 'Too many genes'
    sys.exit(1)
  gene_name = other.columns[0]

  not_amplified = combined_treatment_result[combined_treatment_result[gene_name] < 3]
  amplified = combined_treatment_result[combined_treatment_result[gene_name] >= 3]
  not_amplified_resp = combined_treatment_result[models_treated.name][not_amplified.index].dropna()
  amplified_resp  = combined_treatment_result[models_treated.name][amplified.index].dropna()
  if len(not_amplified_resp) >= 5 and len(amplified_resp) >= 5:
    ttest_res = stats.ttest_ind(not_amplified_resp, amplified_resp, equal_var=False)
    ttest = pd.DataFrame({'t': ttest_res[0], 'p': ttest_res[1]}, index=[gene_name])
    return ttest
  return pd.DataFrame()

def get_cna_deleted_ttest(models_treated, other=None):
  combined_treatment_result = models_treated.T.join(other, how='inner')
  if len(other.columns) != 1:
    print 'Too many genes'
    sys.exit(1)
  gene_name = other.columns[0]

  not_deleted = combined_treatment_result[~combined_treatment_result[gene_name].isin([0,1])]
  deleted = combined_treatment_result[combined_treatment_result[gene_name].isin([0,1])]
  not_deleted_resp  = combined_treatment_result[models_treated.name][not_deleted.index].dropna()
  deleted_resp  = combined_treatment_result[models_treated.name][deleted.index].dropna()
  if len(not_deleted_resp) >= 5 and len(deleted_resp) >= 5:
    ttest_res = stats.ttest_ind(not_deleted_resp, deleted_resp, equal_var=False)
    return pd.DataFrame({'t': ttest_res[0], 'p': ttest_res[1]}, index=[gene_name])
  return pd.DataFrame()

def get_mutation_ttest(models_treated, other=None):
  model_name = models_treated.index[0]
  combined_treatment_result = models_treated.T.join(other, how='inner')
  if len(other.columns) != 1:
    print 'Too many genes'
    sys.exit(1)
  gene_name = other.columns[0]

  not_mutated = combined_treatment_result[combined_treatment_result[gene_name] == 0]
  mutated = combined_treatment_result[combined_treatment_result[gene_name] == 1]
  not_mutated_resp = combined_treatment_result[model_name][not_mutated.index].dropna()
  mutated_resp = combined_treatment_result[model_name][mutated.index].dropna()
  if len(not_mutated_resp) >= 5 and len(mutated_resp) >= 5:
    ttest_res = stats.ttest_ind(not_mutated_resp, mutated_resp, equal_var=False)
    return pd.DataFrame({'t': ttest_res[0], 'p': ttest_res[1]}, index=[gene_name])
  return pd.DataFrame()

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

  annotation['cancer_type'] = annotation['cancer_type'].replace(['GBM', 'LGG'], 'GBMLGG')
  annotation['cancer_type'] = annotation['cancer_type'].replace('COAD/READ', 'COADREAD')
  return annotation

def cna(limitation_row, gene, treatments, cnas_of_type, cnas):
  if gene in cnas.columns:
    pancan_cna_gene_data = cnas[[gene]]
  else:
    print limitation_row.values, 'not in CNA'
    return pd.DataFrame(), pd.DataFrame()

  if gene in cnas_of_type.columns:
    cancer_type_cna_gene_data = cnas_of_type[[gene]]
  else:
    cancer_type_gene_date = pd.DataFrame()

  ct_ttest = pd.DataFrame()
  if limitation_row['Alteration'] == 'amplification':
    ct_ttest = treatments.groupby(level=0).apply(get_cna_amplified_ttest,
                                                 other=cancer_type_cna_gene_data)
    pancan_ttest = treatments.groupby(level=0).apply(get_cna_amplified_ttest,
                                                     other=pancan_cna_gene_data)
  if limitation_row['Alteration'] == 'deletion':
    ct_ttest = treatments.groupby(level=0).apply(get_cna_deleted_ttest,
                                                 other=cancer_type_cna_gene_data)
    pancan_ttest = treatments.groupby(level=0).apply(get_cna_deleted_ttest,
                                                     other=pancan_cna_gene_data)
  return ct_ttest, pancan_ttest

def run(treatments, indir, limiting, outdir):
  annotation = prep_annotation(indir)
  cnas = prep_cnas(indir)
  cnas_with_type = cnas.join(annotation)
  mutations = prep_mutations(indir)
  mutations_with_type = mutations.join(annotation)

  pancan_df = pd.DataFrame()
  for cancer_type, limits in limiting.groupby('Cancer type'):
    cnas_of_type = cnas[cnas_with_type['cancer_type'] == cancer_type]
    mutations_of_type =  mutations[mutations_with_type['cancer_type'] == cancer_type]

    ct_df = pd.DataFrame()
    for i, row in limits.iterrows():
      gene = row['Gene'].split('_')[1][1:]
      print gene
      if row['Alteration'] in ['amplification', 'deletion']:
        ct_ttest, pancan_ttest = cna(row, gene, treatments, cnas_of_type, cnas)
      else:
        if row['Alteration'] == 'mutation':
          if gene in mutations:
            pancan_mut_gene_data = mutations[[gene]]
            pancan_ttest = treatments.groupby(level=0).apply(get_mutation_ttest,
                                                           other=pancan_mut_gene_data)
          else:
            print row.values, 'not in Mutations'

          if gene in mutations_of_type:
            cancer_type_mut_gene_data = mutations_of_type[[gene]]
            ct_ttest = treatments.groupby(level=0).apply(get_mutation_ttest,
                                                         other=cancer_type_mut_gene_data)
        pancan_ttest = pancan_ttest.dropna()
      pancan_ttest['Gene'] = '\'' + gene
      pancan_ttest['Alteration'] = row['Alteration']
      pancan_ttest = pancan_ttest.set_index(['Gene', 'Alteration'], append=True)

      ct_ttest = ct_ttest.dropna()
      ct_ttest['Gene'] = '\'' + gene
      ct_ttest['Alteration'] = row['Alteration']
      ct_ttest = ct_ttest.set_index(['Gene', 'Alteration'], append=True)

      # print row.values
      # print ct_ttest
      # print pancan_ttest
      ct_df = ct_df.append(ct_ttest)
      pancan_df = pancan_df.append(pancan_ttest)

    if len(ct_df.index.names) > 1:
      ct_df.index = ct_df.index.droplevel(level=1)
      ct_df = ct_df[ct_df['p'] < 0.1]
    ct_df.to_csv(os.path.join(outdir, cancer_type + '_cell_lines_limited_ttest.csv'))

    ct_df = pd.DataFrame()
  pancan_df = pancan_df.drop_duplicates()
  pancan_df.index = pancan_df.index.droplevel(level=1)
  pancan_df = pancan_df[pancan_df['p'] < 0.1]
  pancan_df.to_csv(os.path.join(outdir, 'pancan_cell_lines_limited_ttest.csv'))



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

def prep_mutations(indir):
  mutations = pd.read_csv(os.path.join(indir, 'PANCAN_SEQ_BEM.txt'), sep='\t', index_col=0)
  mutations = mutations.T
  mutations.index.rename('COSMIC', inplace=True)
  return mutations


def main():
  indir, limiting, outdir = get_options()

  treatments = prep_treatment(indir)
  limiting = pd.read_csv(limiting)
  run(treatments, indir, limiting, outdir)


if __name__ == "__main__":
  main()

