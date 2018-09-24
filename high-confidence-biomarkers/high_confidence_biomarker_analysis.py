#!/usr/bin/env python
# encoding: utf-8
'''
high_confidence_biomarker_analysis.py

Created by Joan Smith
on 2017-11-04.

Given the collected zscores files created from high_confidence_files script, apply arbitrary criteria and write a file for each one.

Copyright (c) 2018. All rights reserved.
'''
import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob

sys.path.append('../common/')
import utilities as util

TCGA_DIR = '~/Dropbox/genomic-features-survival/gene-only-time-since-surgery/mutation-time-since-surgery-2percent'
CANCER_GENE_FILE = '~/Dropbox/genomic-features-survival/High-confidence biomarkers/cancer_gene_list.csv'

def get_options():
  parser = argparse.ArgumentParser(description='Output univariate cox from pcna25')
  parser.add_argument('-i', action='store', dest='indir')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-s', action='store_true', dest='split_cancer_files')
  namespace = parser.parse_args()

  return (namespace.indir, namespace.output_directory, namespace.split_cancer_files)

def stouffer_sig_and_two_zscore_sig(row):
  if np.absolute(row['Stouffers']) > 3.3:
    row = row.drop('Stouffers')
    if (row > 1.96).sum() >= 2 or (row < -1.96).sum() >= 2:
      return True
  return False

def add_cancer_gene_col(row, cancer_genes=None):
  if row.name == None:
    return None
  if row.name in cancer_genes['cancer_genes'].values:
    row['cancer_gene'] = '***'
  else:
    row['cancer_gene'] = ''
  return row


def apply_criteria(f, criteria, cna_or_mut='cna'):
  df = pd.read_csv(f, index_col=0)
  processing_df = df
  if 'Chromosome' in df.columns:
    processing_df = df.drop(['Chromosome', 'Chromosome Band', 'Location'], axis=1)
  processing_df = processing_df.astype(float)
  criteria_met = processing_df.apply(stouffer_sig_and_two_zscore_sig, axis=1)

  tcga_col = [i for i in df.columns if 'TCGA' in i][0].split('_', 1)[1]
  tcga_glob = os.path.expanduser(os.path.join(TCGA_DIR, tcga_col, tcga_col + '*.zscores.out.csv'))
  tcga_path = glob.glob(tcga_glob)[0]
  tcga_data = pd.read_csv(tcga_path, index_col=0)
  metadata = df[criteria_met].join(tcga_data[['num mutations', 'num patients']], how='left')

  cancer_genes = pd.read_csv(CANCER_GENE_FILE, header=None, names=['cancer_genes'])
  cancer_genes['cancer_genes'] = '\'' + cancer_genes['cancer_genes']
  metadata = metadata.apply(add_cancer_gene_col, cancer_genes=cancer_genes, axis=1)
  return metadata



def main():
  indir, outdir, split_files = get_options()
  files = os.listdir(os.path.join(indir, 'cnas'))
  files = util.remove_extraneous_files(files)

  criteria_list = [stouffer_sig_and_two_zscore_sig]

  for criteria in criteria_list:
    if not split_files:
      outfile_c = open(os.path.join(outdir, criteria.__name__ + '.out.csv'), 'w')
      outfile_m = outfile_c
    for f in files:
      cancer_type = f.split('_')[0]
      if split_files:
        outfile_c = os.path.join(outdir, cancer_type + '_CNA_.criteria_met.out.csv')
        outfile_m = os.path.join(outdir, cancer_type + '_MUT_.criteria_met.out.csv')

      cna_cancer_type_criteria_met = apply_criteria(os.path.join(indir, 'cnas', f),
                                                    criteria, 'cna')
      cna_cancer_type_criteria_met.index = 'cna_' + cna_cancer_type_criteria_met.index
      cna_cancer_type_criteria_met.to_csv(outfile_c, index_label='CNA_'+cancer_type)

      mut_cancer_type_criteria_met = apply_criteria(os.path.join(indir, 'mutations', f),
                                                   criteria, 'mut')
      print f, mut_cancer_type_criteria_met.index
      mut_cancer_type_criteria_met.index = 'mut_' + mut_cancer_type_criteria_met.index
      mut_cancer_type_criteria_met.to_csv(outfile_m, index_label='MUT_'+cancer_type)
  if not split_files:
    outfile.close()


if __name__ == "__main__":
  main()


