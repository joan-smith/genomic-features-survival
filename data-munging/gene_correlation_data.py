#!/usr/bin/env python
# encoding: utf-8
'''
gene_correlation_data.py

Created by Joan Smith
on 2017-6-30

Given a list of genes* from each class of data, collect the original data (with survival information)
for each gene for each patient, for each cancer type, for each class of data

Copyright (c) 2017 . All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np

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

def get_file_class(f):
  return f.split('/')[1]

def collect_data_by_type(cancer_type, data_class_files, gene_list, clinical_dir):
  print cancer_type
  all_classes_data = []
  clinical_file = os.path.join('.', 'clinical', cancer_type + '.clin.merged.txt')
  print clinical_file
  clinical = util.get_clinical_data(clinical_file)
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
      all_classes_data.append(gene_data)

  all_classes_df = pd.concat(all_classes_data, axis=1, verify_integrity=True)
  all_classes_df = all_classes_df.join(clinical, how='left')
  all_classes_df.index.name = cancer_type
  return all_classes_df

def collect_cancer_type_files(cancer_type, indir):
  cancer_type_regex = os.path.join(indir, '*' , '*' + cancer_type + '*')
  cancer_type_files = glob.glob(cancer_type_regex)

  # switch out original cnv data for processed cnv data
  cancer_type_files = [f for f in cancer_type_files if not 'CNV' in f]
  cancer_type_files.extend(glob.glob(os.path.join('.', 'CNV-for-genes', cancer_type + '*')))
  return cancer_type_files


def main():
  indir, clinical_dir, outdir, gene_list_file = get_options()

  clinical_files  = glob.glob(os.path.join(clinical_dir, '*.txt'))
  clinical_files = util.remove_extraneous_files(clinical_files)
  gene_list = pd.read_csv(gene_list_file)

  all_data = []
  for clinical_file in clinical_files:
    cancer_type = util.get_cancer_type(clinical_file)
    cancer_type_files = collect_cancer_type_files(cancer_type, indir)
    cancer_type_data = collect_data_by_type(cancer_type, cancer_type_files, gene_list, clinical_dir)
    all_data.append(cancer_type_data)
  with open(gene_list_file.split('.')[0] + '.output.csv', 'w') as output:
    for cancer_type_data in all_data:
      cancer_type_data.to_csv(output)





if __name__ == "__main__":
  main()

