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
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('-i', action='store', dest='input_directory', required=True)
  parser.add_argument('-l', action='store', dest='limited_query', required=True)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.limited_query, namespace.output_directory)

def get_cna_amplified_ttest(models_treated, other=None):
  if models_treated.name in ['ArmLevelCNScore', 'FocalCNScore']:
    return pd.DataFrame()

  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')

  if len(other.columns) != 1:
    print 'Too many genes'
    sys.exit(1)
  gene_name = other.columns[0]

  not_amplified = combined_treatment_result[combined_treatment_result[gene_name] < 0.3]
  amplified = combined_treatment_result[combined_treatment_result[gene_name] >= 0.3]
  not_amplified_resp  = combined_treatment_result['BestAvgResponse'][not_amplified.index].dropna()
  amplified_resp  = combined_treatment_result['BestAvgResponse'][amplified.index].dropna()
  if len(not_amplified_resp) >= 5 and len(amplified_resp) >= 5:
    if amplified_resp.mean() < 20:
      ttest_res = stats.ttest_ind(not_amplified_resp, amplified_resp, equal_var=False)
      return pd.DataFrame({'t': ttest_res[0], 'p': ttest_res[1],
                        'non mean resp': not_amplified_resp.mean(), 'mean resp': amplified_resp.mean(),
                        'non count': not_amplified_resp.count(), 'count': amplified_resp.count()},
                        index=[gene_name])
  return pd.DataFrame()

def get_cna_deleted_ttest(models_treated, other=None):
  if models_treated.name in ['ArmLevelCNScore', 'FocalCNScore']:
    return pd.DataFrame()

  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')

  if len(other.columns) != 1:
    print 'Too many genes'
    sys.exit(1)
  gene_name = other.columns[0]

  not_deleted = combined_treatment_result[combined_treatment_result[gene_name] > -0.3]
  deleted = combined_treatment_result[combined_treatment_result[gene_name] <= -0.3]
  not_deleted_resp = combined_treatment_result['BestAvgResponse'][not_deleted.index].dropna()
  deleted_resp  = combined_treatment_result['BestAvgResponse'][deleted.index].dropna()
  if len(not_deleted) >= 5 and len(deleted) >= 5:
    if deleted_resp.mean() < 20:
      ttest_res = stats.ttest_ind(not_deleted_resp, deleted_resp, equal_var=False)
      return pd.DataFrame({'t': ttest_res[0], 'p': ttest_res[1],
                        'non mean resp': not_deleted_resp.mean(), 'mean resp': deleted_resp.mean(),
                        'non count': not_deleted_resp.count(), 'count': deleted_resp.count()},
                        index=[gene_name])
  return pd.DataFrame()

def get_mutation_ttest(models_treated, other=None):
  treatment_result = models_treated[['Model', 'BestAvgResponse']]
  treatment_result = treatment_result.set_index('Model')
  combined_treatment_result = treatment_result.join(other, how='inner')
  return pd.DataFrame()

  ttests = {}
  if len(other.columns) != 1:
    print 'Too many genes'
    sys.exit(1)
  gene_name = other.columns[0]

  not_mutated = combined_treatment_result[combined_treatment_result[gene_name] == 0].index
  not_mutated_resp = combined_treatment_result['BestAvgResponse'][not_mutated].dropna()
  mutated = combined_treatment_result[combined_treatment_result[gene_name] == 1].index
  mutated_resp = combined_treatment_result['BestAvgResponse'][mutated].dropna()
  if len(not_mutated_resp) >= 5 and len(mutated_resp) >= 5:
    ttest_res = stats.ttest_ind(not_mutated_resp, mutated_resp, equal_var=False)
    return pd.DataFrame({'t': ttest_res[0], 'p': ttest_res[1],
                        'non mean resp': not_mutated_resp.mean(), 'mean resp': mutated_resp.mean(),
                        'non count': not_mutated_resp.count(), 'count': mutated_resp.count()},
                        index=[gene_name])
  return pd.DataFrame()

def get_copy_number(gene):
  copy_number = gene.str.split(',', expand=True)[0]
  copy_number = copy_number.astype(int)
  copy_number = copy_number.replace(-1, np.nan)
  return copy_number

def prep_cnas(xls):
  cna = xls.parse('copy number', index_col=0).T
  cna = cna.divide(2)
  cna = cna.apply(np.log2)
  return cna

def cna(limitation_row, gene, treatments, cnas):
  if gene in cnas.columns:
    pancan_cna_gene_data = cnas[[gene]]
  else:
    return pd.DataFrame()

  if limitation_row['Alteration'] == 'amplification':
    pancan_ttest = treatments.groupby(['Treatment']).apply(get_cna_amplified_ttest,
                                                     other=pancan_cna_gene_data)
  if limitation_row['Alteration'] == 'deletion':
    pancan_ttest = treatments.groupby(['Treatment']).apply(get_cna_deleted_ttest,
                                                     other=pancan_cna_gene_data)
  return pancan_ttest

def run(xls, limits, outdir):
  cnas = prep_cnas(xls)
  mutations = prep_mutations(xls)
  treatments = prep_treatment(xls)

  pancan_df = pd.DataFrame()
  for i, row in limits.iterrows():
    gene = row['Gene'].split('_')[1][1:]
    print gene
    if row['Alteration'] in ['amplification', 'deletion']:
      pancan_ttest = cna(row, gene, treatments, cnas)
    else:
      if row['Alteration'] == 'mutation':
        if gene in mutations:
          pancan_mut_gene_data = mutations[[gene]]
          pancan_ttest = treatments.groupby(['Treatment']).apply(get_mutation_ttest,
                                                         other=pancan_mut_gene_data)

    pancan_ttest = pancan_ttest.dropna()
    pancan_ttest['Gene'] = '\'' + gene
    pancan_ttest['Alteration'] = row['Alteration']
    pancan_ttest = pancan_ttest.set_index(['Gene', 'Alteration'], append=True)

    print row.values
    print pancan_ttest
    pancan_df = pancan_df.append(pancan_ttest)

  pancan_df = pancan_df.drop_duplicates()
  pancan_df.index = pancan_df.index.droplevel(level=1)
  pancan_df = pancan_df[pancan_df['p'] < 0.1]
  pancan_df.to_csv(os.path.join(outdir, 'pancan_xenografts_limited_ttest.csv'))



def clean_treatment(t):
  t = t.drop('COSMIC', axis=1)
  return t.mean()


def prep_xls(infile):
  xls = pd.ExcelFile(infile)
  return xls

def prep_treatment(xls):
  results = xls.parse('PCT curve metrics')
  return results

def prep_mutations(xls):
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



def main():
  indir, limiting, outdir = get_options()

  xls = prep_xls(os.path.join(indir, 'nm.3954-S2.xlsx'))
  limiting = pd.read_csv(limiting)
  run(xls, limiting, outdir)


if __name__ == "__main__":
  main()

