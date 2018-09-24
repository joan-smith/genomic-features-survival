#!/usr/bin/env python
# encoding: utf-8
'''
drug_sensitivity_cell_lines.py

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
from IPython import embed

def get_options():
  parser = argparse.ArgumentParser(description='Calculate drug sensitivity correlations for cell lines')
  parser.add_argument('-i', action='store', dest='input_dir', required=True)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_dir, namespace.output_directory)

def get_correlations(treatment, other=None):
  treatment = treatment.T
  treatment.columns = ['Response']
  combined_treatment_result = other.join(treatment, how='inner')
  corrs = combined_treatment_result.corrwith(combined_treatment_result['Response'])
  return corrs


def prep_annotation(indir):
  annotation_xls = pd.ExcelFile(os.path.join(indir, 'TableS1E - annotations.xlsx'))
  annotation = annotation_xls.parse('TableS1E-CellLines', parse_cols='C,K', header=1)
  annotation.columns = ['COSMIC', 'cancer_type']
  annotation = annotation.dropna(how='all')
  # COSMIC is being parsed as a float, but it needs to be a string.
  annotation['COSMIC'] = annotation['COSMIC'].astype(int).astype(str)
  annotation = annotation.set_index('COSMIC')
  return annotation

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

def run_mutations_analysis(annotation, treatments, mutations):
  pancan_correlations = treatments.groupby(level=0).apply(get_correlations,
                                  other=mutaions)
  pancan_correlations.T.to_csv('pancan_correlations.csv')

  cancer_type_groups = mutations.join(annotation).groupby('cancer_type')
  for cancer_type, mutations_of_type in cancer_type_groups:
    cancer_type = cancer_type.replace('/', '')
    mutations_of_type = mutations_of_type.drop('cancer_type', axis=1)
    print cancer_type, mutations_of_type.shape[0]
    if mutations_of_type.shape[0] >= 10:
      corr = treatments.groupby(level=0).apply(get_correlations,
                                               other=mutations_of_type)
      print corr
      corr.T.to_csv(cancer_type + '_correlations.csv')


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


def run_cna(treatments, annotations, outdir):
  cnas = prep_cnas(indir)
  print cnas

  cancer_type_groups = cnas.join(annotation).groupby('cancer_type')
  for cancer_type, cnas_of_type in cancer_type_groups:
    cancer_type = cancer_type.replace('/', '')
    cnas_of_type = cnas_of_type.drop('cancer_type', axis=1)
    print cancer_type, cnas_of_type.shape[0]
    if cnas_of_type.shape[0] >= 10:
      corr = treatments.groupby(level=0).apply(get_correlations,
                                               other=cnas_of_type)
      corr.T.to_csv(os.path.join(outdir, cancer_type + '_correlations.csv'))


def do_work(indir, outdir):
  print indir
  annotation = prep_annotation(indir)
  treatments = prep_treatment(indir)
  # mutations = prep_mutations(indir)
  # run_mutations_analysis(annotation, treatments, mutations)




def main():
  indir, outdir = get_options()
  do_work(indir, outdir)

  treatments = prep_treatment(indir)
  cnas = prep_cnas(indir)
  cnas.T.to_csv('raw_cnas.csv')

  pancan_correlations = treatments.groupby(level=0).apply(get_correlations,
                                  other=cnas)
  pancan_correlations.T.to_csv('pancan_cna_correlations.csv')




if __name__ == "__main__":
  main()

