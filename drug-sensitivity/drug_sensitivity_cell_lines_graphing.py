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


def run(treatments, indir, limiting, outdir):
  annotation = prep_annotation(indir)
  cnas = prep_cnas(indir)
  cnas_with_type = cnas.join(annotation)
  mutations = prep_mutations(indir)
  mutations_with_type = mutations.join(annotation)

  pancan_df = pd.DataFrame()
  for i, row in limiting.iterrows():
    row = row.str.strip()
    file_name = '_'.join(row) + '.csv'
    print row
    drug_treatments = treatments.loc[[row['Drug']]].T
    if row['Aberration'] == 'MUT':
      if row['Cancer type'] in ['Pan-cancer']:
        mutations_of_type = mutations
      else:
        mutations_of_type =  mutations[mutations_with_type['cancer_type'] == row['Cancer type']]
      gene_data =  mutations_of_type[row['Gene']]
    if row['Aberration'] == 'CNA':
      if row['Cancer type'] in ['Pan-cancer']:
        cnas_of_type = cnas
      else:
        cnas_of_type = cnas[cnas_with_type['cancer_type'] == row['Cancer type']]
      gene_data = cnas_of_type[row['Gene']]
    gene_drug_treatments = drug_treatments.join(gene_data, how='inner')
    gene_drug_treatments.to_csv(os.path.join(outdir, file_name))

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

