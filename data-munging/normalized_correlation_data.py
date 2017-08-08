#!/usr/bin/env python
# encoding: utf-8
'''
normalized_correlation_data.py

Created by Joan Smith
on 2017-7-21

Given a list of genes* from each class of data, collect the normalized data
for each gene for each patient, for each cancer type, for each class of data.
Output both the normalized data and the pearson correlation matrix

Copyright (c) 2017 . All rights reserved.
'''

import argparse
import sys
import os
import glob
import itertools
import pandas as pd
import numpy as np
import scipy.stats as stats
import pdb

sys.path.append('../common/')
import utilities as util
import analysis

SEP_BY_CLASS = {
  'CNV-for-genes': ',',
  'RNASeq-data': '\t',
  'RPPA-data': ',',
  'methylation-mean': '\t',
  'microRNA': ',',
  'mutation-data': '\t',
}

def get_options():
  parser = argparse.ArgumentParser(description='Produce intermediate file for analysis.')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-l', action='store', dest='gene_list_file')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory,
    namespace.output_directory, namespace.gene_list_file)

def subtract_mean_norm(df):
  means = df.mean(axis=0)
  print 'means'
  print means
  normed = df.sub(means, axis=1)
  return normed

def rnaseq_norm(df):
  df_log2 = df.apply(np.log2)
  clipped = np.clip(df_log2, 0, np.inf)
  means = clipped.mean(axis=0)
  print 'means', means
  normed = clipped.sub(means, axis=1)
  return normed

def microrna_norm(df):
  clipped = np.clip(df, 0, np.inf)
  means = clipped.mean(axis=0)
  print 'means', means
  normed = clipped.sub(means, axis=1)
  return normed


def no_norm(df):
  return df


NORMALIZATIONS = {
  'CNV-for-genes': subtract_mean_norm,
  'RNASeq-data': rnaseq_norm,
  'RPPA-data': subtract_mean_norm,
  'methylation-mean': subtract_mean_norm,
  'microRNA': microrna_norm,
  'mutation-data': no_norm,
}


def get_file_class(f):
  return f.split('/')[1]


def collect_data_by_type(cancer_type, data_class_files, gene_list):
  print cancer_type
  all_classes_data = []
  for f in data_class_files:
    file_class = get_file_class(f)
    print 'File Class:', file_class
    df = pd.read_csv(f, sep=SEP_BY_CLASS[file_class], index_col=0, low_memory=False)
    if file_class in gene_list:
      requested_genes = gene_list[file_class]
      gene_data = df.loc[requested_genes.dropna()]
      gene_data = gene_data.transpose().reset_index()
      if file_class  == 'CNV-for-genes':
        gene_data = gene_data.set_index('index')
        gene_data = gene_data.drop(['Chromosome', 'Location'])
      else:
        gene_data = util.maybe_clear_non_01s(gene_data, 'index', cancer_type)
        gene_data = util.add_identifier_column(gene_data, 'index')
        gene_data = gene_data.drop('index', axis=1)
        gene_data = gene_data.set_index('identifier')
      gene_name = lambda x: '\'' + x if not '\'' in x else x
      gene_data.rename(columns=lambda x: file_class + ':' + gene_name(x), inplace=True)
      gene_data = gene_data.astype(float)
      normed_gene_data = NORMALIZATIONS[file_class](gene_data)
      all_classes_data.append(normed_gene_data)

  all_classes_df = pd.concat(all_classes_data, axis=1, verify_integrity=True)
  all_classes_df.index.name = cancer_type
  return all_classes_df

def collect_cancer_type_files(cancer_type, indir):
  cancer_type_regex = os.path.join(indir, '*' , '*' + cancer_type + '*')
  cancer_type_files = glob.glob(cancer_type_regex)

  # switch out original cnv data for processed cnv data
  cancer_type_files = [f for f in cancer_type_files if not 'CNV' in f]
  cancer_type_files.extend(glob.glob(os.path.join('.', 'CNV-for-genes', cancer_type + '*')))
  return cancer_type_files

def run_pearsons(all_data_df):
  combinations = itertools.combinations(all_data_df.columns, 2)
  print combinations

  correlations = pd.DataFrame(index=all_data_df.columns)
  for pair in combinations:
    row, col = pair
    temp_df = pd.DataFrame([all_data_df[row], all_data_df[col]])
    temp_df = temp_df.dropna(how='any', axis=1)
    c, p = stats.pearsonr(temp_df.loc[row].values, temp_df.loc[col].values)
    correlations.loc[row, col] = c
    correlations.loc[col, row] = c
    print correlations.loc[row, col]
  return correlations

def main():
  indir, clinical_dir, outdir, gene_list_file = get_options()
  input_name = os.path.basename(gene_list_file).split('.')[0]
  normed_output_file = os.path.join(outdir,  input_name + '.normed_output.csv')
  correlation_output_file = os.path.join(outdir,  input_name + '.correlation_output.csv')

  clinical_files  = glob.glob(os.path.join(clinical_dir, '*.txt'))
  clinical_files = util.remove_extraneous_files(clinical_files)
  gene_list = pd.read_csv(gene_list_file)

  all_data = []
  for clinical_file in clinical_files:
    cancer_type = util.get_cancer_type(clinical_file)
    cancer_type_files = collect_cancer_type_files(cancer_type, indir)
    cancer_type_data = collect_data_by_type(cancer_type, cancer_type_files, gene_list)
    all_data.append(cancer_type_data)

  with open(normed_output_file, 'w') as output:
    for cancer_type_data in all_data:
      cancer_type_data.to_csv(output)

  all_data_df = pd.concat(all_data, verify_integrity=True)
  correlations = run_pearsons(all_data_df)
  correlations.to_csv(correlation_output_file)




if __name__ == "__main__":
  main()

